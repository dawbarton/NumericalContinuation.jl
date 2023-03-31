mutable struct Chart{T}
    const u::Vector{T}
    const ts::Vector{T}
    s::Int
    h::T
end

Base.copy(chart::Chart) = Chart(copy(chart.u), copy(chart.ts), chart.s, chart.h)

function copyto_chart!(dest::Chart, src::Chart)
    copyto!(dest.u, src.u)
    copyto!(dest.ts, src.ts)
    dest.s = src.s
    dest.h = src.h
    return dest
end

struct PseudoArclengthFunc{T}
    chart::Chart{T}
    tmp::DiffCache{Vector{T}, Vector{T}}
end

function PseudoArclengthFunc(chart::Chart)
    return PseudoArclengthFunc(chart, DiffCache(similar(chart.u)))
end

function (paf::PseudoArclengthFunc)(u)
    tmp = get_tmp(paf.tmp, u)  # avoid allocations
    tmp .= u .- paf.chart.u
    return dot(tmp, paf.chart.ts)
end

struct PseudoArclength{T, S}
    h0::T
    h_min::T
    h_max::T
    h_grow::T
    h_shrink::T
    alpha_max::T
    maxiters::Int
    abstol::T
    nl_solver::S
end

function PseudoArclength(;
                         h0::Real = 1 // 1000,
                         h_min::Real = 1 // 1_000_000,
                         h_max::Real = 1,
                         h_grow::Real = 6 // 5,
                         h_shrink::Real = 1 // 2,
                         alpha_max::Real = 99 // 100,
                         maxiters::Integer = 1000,
                         abstol::Real = 1 // 1_000_000,
                         nl_solver = NewtonRaphson())
    _h0, _h_min, _h_max, _h_grow, _h_shrink, _alpha_max, _abstol = promote(h0, h_min, h_max,
                                                                           h_grow, h_shrink,
                                                                           alpha_max,
                                                                           abstol)
    PseudoArclength{typeof(_h0), typeof(nl_solver)}(_h0, _h_min, _h_max, _h_grow, _h_shrink,
                                                    _alpha_max, maxiters, _abstol,
                                                    nl_solver)
end

struct PseudoArclengthCache{T, F, S}
    prob::ContinuationProblem
    f::F
    nl_solver::S
    h0::T
    h_min::T
    h_max::T
    h_grow::T
    h_shrink::T
    alpha_max::T
    maxsteps::Int
    chart_base::Chart{T}  # Corrected
    chart_next::Chart{T}  # May not be corrected; TODO: is chart_next_status required?
    atlas::Vector{Chart{T}}
    step::Base.RefValue{Int}
    force_stop::Base.RefValue{Bool}
end

function SciMLBase.init(prob::ContinuationProblem, alg::PseudoArclength,
                        args...; pars = nothing, maxsteps::Integer = 100,
                        dir::Integer = 1, kwargs...)
    # Copy the problem wrapper to avoid mutating the original
    pa_prob = ContinuationProblem(prob.zero_function!, copy(prob.monitor_function),
                                  copy(prob.monitor_function_name), prob.sub_problem,
                                  prob.sub_problem_name)
    # Convert all the parameters to the same type
    T = get_eltype(pa_prob)
    h0 = convert(T, alg.h0)
    h_min = convert(T, alg.h_min)
    h_max = convert(T, alg.h_max)
    h_grow = convert(T, alg.h_grow)
    h_shrink = convert(T, alg.h_shrink)
    alpha_max = convert(T, alg.alpha_max)
    # Check the direction is +/- 1
    if dir != 1 && dir != -1
        throw(ArgumentError("dir must be 1 or -1"))
    end
    # Add the pseudo-arclength condition
    n = get_eqns(pa_prob) + 1
    chart_base = Chart(zeros(T, n), zeros(T, n), dir, h0)
    chart_next = Chart(zeros(T, n), zeros(T, n), dir, h0)
    add_monitor_function!(pa_prob, :pseudoarclength,
                          monitor_function(PseudoArclengthFunc(chart_next); active = false,
                                           value = 0))
    # Construct the function wrapper
    f = ContinuationWrapper(pa_prob, pars)
    (u0, p0) = get_initial(f)
    if length(u0) != n
        throw(DimensionMismatch("Equations and unknowns must match: currently $n equations and $(length(u0)) unknowns"))
    end
    chart_next.u .= u0
    chart_next.ts[end] = one(T)
    # Construct the nonlinear solver
    nl_solver = init(NonlinearProblem{true}(f, u0, p0), alg.nl_solver)
    return PseudoArclengthCache{T, typeof(f), typeof(nl_solver)}(pa_prob, f, nl_solver,
                                                                 h0, h_min, h_max, h_grow,
                                                                 h_shrink, alpha_max,
                                                                 maxsteps, chart_base,
                                                                 chart_next, Chart{T}[],
                                                                 Ref(0), Ref(false))
end

function SciMLBase.solve!(cache::PseudoArclengthCache)
    initial_correct!(cache)
    cache.force_stop[] || save_chart!(cache)
    cache.step[] += 1
    while !cache.force_stop[] && (cache.step[] < cache.maxsteps)
        predict!(cache)
        cache.force_stop[] || correct!(cache)
        cache.force_stop[] || save_chart!(cache)
        cache.step[] += 1
    end
    return cache.atlas
end

function initial_correct!(cache::PseudoArclengthCache)
    (; nl_solver, chart_base, chart_next) = cache
    # Initial correction
    reinit!(nl_solver, chart_next.u)
    sol = solve!(nl_solver)
    if sol.retcode != ReturnCode.Success
        @debug "Failed to converge on initial correction"
        cache.force_stop[] = true
    else
        chart_next.u .= sol.u
        # Calculate tangent space
        tangent_space!(chart_next.ts, nl_solver)
        if chart_next.ts[end] < 0
            chart_next.ts .= .-chart_next.ts
        end
        # Save successful chart as the new base
        copyto_chart!(chart_base, chart_next)
    end
    return
end

function predict!(cache::PseudoArclengthCache)
    (; chart_base, chart_next, h_min) = cache
    if chart_next.h < h_min
        cache.force_stop = true
        return
    end
    # Linear predicter
    chart_next.u .= chart_base.u .+ chart_next.s .* chart_next.h .* chart_base.ts
    return
end

function tangent_space!(ts, nl_solver)
    NonlinearSolve.jacobian!(nl_solver.J, nl_solver)
    for i in axes(nl_solver.J, 2)
        nl_solver.J[end, i] = 0  # zero out the pseudo-arclength condition
    end
    null = qr!(nl_solver.J')
    ts .= null.Q[:, end]
    return ts
end

function correct!(cache::PseudoArclengthCache)
    (; nl_solver, chart_base, chart_next, alpha_max, h_grow, h_shrink, h_max) = cache
    # Correct
    reinit!(nl_solver, chart_next.u)
    sol = solve!(nl_solver)
    # Check for convergence
    if sol.retcode != ReturnCode.Success
        @debug "Failed to converge; retrying with smaller stepsize" i
        # Failed to converge, try again with smaller stepsize
        chart_next.h *= h_shrink
    else
        chart_next.u .= sol.u
        # Calculate tangent space
        tangent_space!(chart_next.ts, nl_solver)
        # Check that the angle between consecutive tangent spaces is sufficiently small
        alpha = dot(chart_base.ts, chart_next.ts)
        if alpha < 0
            # Make sure that the handedness of the nullspace matches the previous one
            chart_next.ts .= .-chart_next.ts
            alpha = -alpha
        end
        if alpha < alpha_max
            # Angle was too high - shrink stepsize and retry
            @debug "Angle between consecutive tangent spaces too large on step; retrying with smaller stepsize" i angle=acosd(alpha)
            chart_next.h *= h_shrink
        else
            # Save successful chart as the new base
            copyto_chart!(chart_base, chart_next)
            # If angle is small (alpha ≈ 1), increase stepsize up to h_max; otherwise leave untouched
            chart_next.h = alpha > (alpha_max + 1) / 2 ?
                           min(chart_next.h * h_grow, h_max) : chart_next.h
        end
    end
    return
end

function save_chart!(cache::PseudoArclengthCache)
    (; chart_base, atlas) = cache
    # Update atlas
    push!(atlas, copy(chart_base))
    return
end
