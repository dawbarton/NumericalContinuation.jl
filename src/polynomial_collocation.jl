struct PolynomialCollocation{T, F}
    f!::F
    u0::Any
    vars::Int
    eqns::Int
    ndim::Int
    nmesh::Int
    ncoll::Int
    repr_poly::OrthogonalPolynomial{T}
    coll_poly::OrthogonalPolynomial{T}
    repr_pts::Vector{T}
    coll_pts::Vector{T}
    In::Matrix{T}
    Dt::Matrix{T}
    InDt::Matrix{T}
    phase::Bool
    u_tmp::DiffCache{Vector{T}, Vector{T}}  # Temporary storage for u
    h::Vector{T}  # Mesh spacing (chart data)
    phase_du::Matrix{T} # Time derivative of u (chart data)
end

function PolynomialCollocation(f, u, tspan, p = (); phase::Bool = true,
                               ncoll::Integer = 4,
                               nmesh::Union{Missing, Integer} = missing,
                               repr = equispaced, coll = legendre)
    length(tspan) == 2 || throw(ArgumentError("tspan must contain two values"))
    if u isa AbstractArray
        ndims(u) == 2 || throw(ArgumentError("u must be a matrix"))
        ndim = size(u, 1)
        mod(size(u, 2) - 1, ncoll) == 0 ||
            throw(ArgumentError("size(u, 2) must be nmesh * ncoll + 1"))
        ismissing(nmesh) || @warn "nmesh is ignored when u is a matrix"
        nmesh = (size(u, 2) - 1) ÷ ncoll
        T = eltype(u)
    else
        u0 = u(tspan[begin])
        ndim = length(u0)
        ismissing(nmesh) && throw(ArgumentError("nmesh must be specified"))
        T = eltype(u0)
    end
    _f = IIPWrapper(f, 3)  # Ensure that the function has an in-place form
    _vars = (nmesh * ncoll + 1) * ndim + length(p) + 2
    _eqns = nmesh * ncoll * ndim + phase  # doesn't include the boundary conditions

    repr_poly = repr(T, ncoll)
    if !is_closed_interval(repr_poly)
        throw(ArgumentError("The representation polynomial must be closed (i.e., the endpoints are included as a node)"))
        # This could be relaxed by an additional interpolation step but doesn't seem worth it at the moment
    end
    repr_pts = (nodes(repr_poly) .+ 1) ./ 2
    coll_poly = coll(T, ncoll - 1)
    coll_pts = (nodes(coll_poly) .+ 1) ./ 2
    In = collect(transpose(interpolation_matrix(repr_poly, nodes(coll_poly))))
    Dt = collect(transpose(differentiation_matrix(coll_poly) .* 2))
    InDt = In * Dt
    h = fill(one(T) / nmesh, nmesh)

    if u isa AbstractArray
        uu = collect(u)::Matrix{T}
        _u = vec(uu)
    else
        uu = Matrix{T}(undef, ndim, nmesh * ncoll + 1)
        t = range(tspan[begin], tspan[end], length = nmesh + 1)
        period = tspan[end] - tspan[begin]
        mesh_t = period / nmesh
        for i in Base.OneTo(nmesh)
            for j in Base.OneTo(ncoll)
                @show _t = t[i] + repr_pts[j] * mesh_t
                uu[:, (i - 1) * ncoll + j + 1] = u(_t)
            end
        end
        uu[:, end] = u(tspan[end])
        _u = vec(uu)
    end

    if phase
        phase_du = Matrix{T}(undef, ndim, nmesh * ncoll)
        i0 = 1
        i1 = ncoll
        for i in Base.OneTo(nmesh)
            @views mul!(phase_du[:, i0:i1], uu[:, i0:(i1 + 1)], InDt)
            i0 += ncoll
            i1 += ncoll
        end
        phase_du ./= norm(phase_du)
    else
        phase_du = zeros(T, ndim, nmesh * ncoll)
    end

    u_tmp = DiffCache(Vector{T}(undef, ndim))

    return PolynomialCollocation{T, typeof(_f)}(_f,
                                                (u = _u, p = p,
                                                 tspan = [tspan[begin], tspan[end]]),
                                                _vars, _eqns, ndim, nmesh, ncoll, repr_poly,
                                                coll_poly, repr_pts, coll_pts, In, Dt, InDt,
                                                phase, u_tmp, h, phase_du)
end

function (coll::PolynomialCollocation{S})(res, uu, data; kwargs...) where {S}
    u = uu.zero
    U = reshape(u.u, (coll.ndim, coll.nmesh * coll.ncoll + 1))
    u_tmp = get_tmp(coll.u_tmp, u.u)
    res_mat = reshape(@view(res[1:(coll.ndim * coll.nmesh * coll.ncoll)]),
                      (coll.ndim, coll.nmesh * coll.ncoll))
    T = u.tspan[end] - u.tspan[begin]
    t0 = u.tspan[begin]
    i0 = 1
    i1 = coll.ncoll
    phase = zero(S)
    for i in Base.OneTo(coll.nmesh)
        for j in Base.OneTo(coll.ncoll)
            # interpolate to find u at the collocation point
            @views mul!(u_tmp, U[:, i0:(i1 + 1)], coll.In[:, j])
            # evaluate the function at the collocation point (multiplied by the period in the next step)
            coll.f!(@view(res_mat[:, i0 + j - 1]), u_tmp, u.p, t0 + coll.repr_pts[j] * T)
            # compute the time derivative of u at the collocation point and subtract from the residual
            @views mul!(res_mat[:, i0 + j - 1], U[:, i0:(i1 + 1)], coll.InDt[:, j], -one(S),
                        T * coll.h[j])
            for k in Base.OneTo(coll.ndim)
                phase += coll.phase_du[k, i0 + j - 1] * u_tmp[k]
            end
        end
        i0 += coll.ncoll
        i1 += coll.ncoll
        t0 += coll.h[i] * T
    end
    if coll.phase
        res[end] = phase
    end
    return res
end

"""
    polynomial_collocation(f, u, tspan, [p];
                           [phase = true], [t0 = 0], [ncol = 4], [nmesh],
                           [repr = equispaced], [coll = legendre])

Implement an orthogonal polynomial-based collocation scheme for discretising a periodic
orbit. The vector field `f` is assumed to be of the same form as used in the SciML/DiffEq
ecosystem, i.e., `f(u, p, t)`, where `u` is the state vector, `p` are the (continuation)
parameters, and `t` is time.

The time interval is discretised into `nmesh` intervals, each of which is represented by a
polynomial of order `ncoll`.

The state `u` is either a matrix or a function `u(t)`.

If `u` is a matrix, the columns are the state vectors at different times. The times are
assumed to equispaced between `tspan[1]` and `tspan[2]`. The number of mesh intervals is
inferred from the number of columns (i.e., size(u, 2) == nmesh * ncoll + 1).

If `u` is a function, it is called with the times at the representation points to return the
state vector at those times. In this case `nmesh` must be specified.

Note that all parameters in `p` are added as (initially inactive) continuation parameters.
If your problem has many parameters that are not likely to be used for continuation, it can
be better to separate them out using a callable `struct` (or a closure) to encapsulate that
data.

By default, the start time `t0` is fixed as zero. This can be changed to be any numerical
value consistent with the problem, or removed entirely by setting `t0 = nothing`.

The representation points are by default equispaced

# Example

Construct an initial solution using a simulation using OrdinaryDiffEq.

```
using OrdinaryDiffEq

function hopf!(res, u, p, t)
    ss = u[1]^2 + u[2]^2
    res[1] = p[1] * u[1] - u[2] + p[2] * u[1] * ss
    res[2] = u[1] + p[1] * u[2] + p[2] * u[2] * ss
    return res
end

# Starting point on the limit cycle with period 2π
u0 = [1.0, 0.0]
tspan = (0.0, 2π)
p = [1.0, -1.0]
odeprob = ODEProblem(hopf!, u0, tspan, p)
sol = solve(odeprob, Tsit5())

# Collocation problem
prob = polynomial_collocation(hopf!, sol, tspan, p)
```
"""
function polynomial_collocation(f, u, tspan, p = (); t0 = 0, kwargs...)
    coll_prob = PolynomialCollocation(f, u, tspan, p; kwargs...)
    prob = ContinuationProblem(coll_prob)
    # Add continuation parameters (all inactive to start with)
    add_parameters!(prob, :p, keys(p))
    # Fix start time
    if t0 !== nothing
        add_parameter!(prob, :t0, @optic(_.tspan[begin]); value = t0)
    end
    return prob
end

struct LimitCycle
    u0::Any
    vars::Int
    eqns::Int
end

LimitCycle(coll::PolynomialCollocation) = LimitCycle((;), 0, coll.ndim)

function (lc::LimitCycle)(res, u, data; kwargs...)
    v = u.coll.zero.u
    res .= v[1:lc.eqns] .- v[end-lc.eqns+1:end]
    return res
end

function limit_cycle(f, u, tspan, p = (); kwargs...)
    coll = polynomial_collocation(f, u, tspan, p; kwargs...)
    lc = LimitCycle(coll.zero_function!)
    return ContinuationProblem(lc; subprob=[:coll => coll])
end

@testitem "Polynomial collocation" begin
    using LinearAlgebra: norm

    function hopf!(res, u, p, t)
        ss = u[1]^2 + u[2]^2
        res[1] = p[1] * u[1] - u[2] + p[2] * u[1] * ss
        res[2] = u[1] + p[1] * u[2] + p[2] * u[2] * ss
        return res
    end

    # Define initial solution
    p0 = [1.0, -1.0]
    t = range(0, 2π, length = 101)
    u0 = [sqrt(p0[1]) .* sin.(t) -sqrt(p0[1]) .* cos.(t)]'
    prob = NumericalContinuation.polynomial_collocation(hopf!, u0, (0, 2π), p0)

    # Problem set up
    (_u, data) = NumericalContinuation.get_initial(prob)
    _active = NumericalContinuation.get_initial_active(prob)
    _monitor = NumericalContinuation.get_initial_monitor(prob, _u, data)
    res_layout = NumericalContinuation.get_initial_residual_layout(prob)
    chart = nothing

    # Buffers
    active = ComponentVector{Bool}(_active)
    u = ComponentVector(ComponentVector{Float64}(_u);
                        monitor = zeros(Float64, count(active)))
    monitor = ComponentVector{Float64}(_monitor)
    res = ComponentVector{Float64}(res_layout)
    @test length(res) == length(u) - size(u0, 1)  # No boundary conditions to complete the problem

    # Optimised code
    func = NumericalContinuation.ContinuationFunction(prob)

    # Evaluate
    NumericalContinuation.eval_zero_function!(res, func, u, data, active, monitor)
    @test norm(res) < 2e-4
end
