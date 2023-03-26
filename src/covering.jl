# Need to shift from IIP to OOP

struct PseudoArclength{U, F}
    f::F
    predicted::U
    tangent::U
end

(pa::PseudoArclength)(u, p) = [pa.f(u, p); dot(u - pa.predicted, pa.tangent)]

function (pa::PseudoArclength)(res, u, p)
    pa.f(@view(res[1:(end - 1)]), u, p)
    res[end] = dot(u - pa.predicted, pa.tangent)
    return res
end

struct Chart{T, U}
    u::U
    ts::U
    s::Int
    h::T
end

"""
    continuation(f, u0, active; [p0], [bounds], [solver])

Continue the zero function `f` in the variable `active`, which is an index into the state
vector `u0`.

The zero function should have a dimensionality deficit of 1. I.e., it should have one more
input than outputs.

`u0` must be a mutable container.
"""
function continuation(f, u0, active; bounds = (), p0 = nothing, solver = TrustRegion(),
                      max_iters = 100, s = 1, h0 = 1e-3, h_min = 1e-8, h_max = 1,
                      h_grow = 1.2, h_shrink = 0.5, alpha_max = 0.99)
    T = eltype(u0)
    if (s != 1) && (s != -1)
        throw(ArgumentError("Initial direction s must be +1 or -1"))
    end
    if !isempty(bounds) && (length(bounds) != 2)
        throw(ArgumentError("bounds must have two elements - an upper and lower value for the active variable"))
    end
    return _continuation(f, u0, active; bounds = bounds, p0 = p0, solver = solver,
                         max_iters = max_iters, s = s, h0 = T(h0), h_min = T(h_min),
                         h_max = T(h_max), h_grow = T(h_grow), h_shrink = T(h_shrink),
                         alpha_max = T(alpha_max))
end

function _continuation(f, u0, active; bounds, p0, solver, max_iters, s, h0, h_min, h_max,
                       h_grow, h_shrink, alpha_max)
    function ts_from_sol(u)
        J = jacobian(f_p0, u, cfg)
        TS = nullspace(J)
        if size(TS, 2) != 1
            throw(ErrorException("Tangent space must be exactly one dimension"))
        end
        return vec(TS)
    end
    # Variables for the pseudo-arclength condition
    predicted = copy(u0)
    tangent = zero(u0)
    tangent[active] = 1
    pa = PseudoArclength(f, predicted, tangent)
    # Create the problem structure for the nonlinear solver
    prob = NonlinearProblem(pa, u0, p0)
    # Correct the initial solution
    sol = solve(prob, solver)
    # Check for convergence
    if sol.retcode != ReturnCode.Success
        throw(ErrorException("Failed to converge on initial solution"))
    end
    # ForwardDiff config
    f_p0 = Base.Fix2(f, p0)
    cfg = JacobianConfig(f_p0, u0)
    # Construct initial chart
    chart = Chart(sol.u, ts_from_sol(sol.u), s, h0)
    atlas = [chart]
    # Iterate
    for i in Base.OneTo(max_iters)
        if chart.h < h_min
            @warn "Minimum stepsize reached"
            break
        end
        # Update pseudo-arclength condition with predicted and tangent space
        pa.predicted .= chart.u .+ chart.s .* chart.h .* chart.ts
        pa.tangent .= chart.ts
        # Correct
        sol = solve(prob, solver; u0 = pa.predicted)
        # Check for convergence
        if sol.retcode != ReturnCode.Success
            # Failed to converge, try again with smaller stepsize
            chart = Chart(chart.u, chart.ts, chart.s, chart.h * h_shrink)
        else
            # Calculate tangent space
            ts = ts_from_sol(sol.u)
            # Check that the angle between consecutive tangent spaces is sufficiently small
            alpha = dot(chart.ts, ts)
            if alpha < 0
                # Make sure that the handedness of the nullspace matches the previous one
                ts .= .-ts
            end
            if alpha < alpha_max
                # Angle was too high - shrink stepsize and retry
                chart = Chart(chart.u, chart.ts, chart.s, chart.h * h_shrink)
            else
                # If angle is small (alpha ≈ 1), increase stepsize up to h_max; otherwise leave untouched
                h = alpha > (alpha_max + 1) / 2 ? min(chart.h * h_grow, h_max) : chart.h
                chart = Chart(sol.u, ts, s, h)
                push!(atlas, chart)
            end
        end
    end
    return atlas
end

function continuation(prob::ContinuationProblem, pars = nothing; kwargs...)
    cw = ContinuationWrapper(prob, pars)
    (u0, p0) = get_initial(cw)

end
