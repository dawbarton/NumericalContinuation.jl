"""
    fourier_diff([T=Float64,] N; order=1)

Create a Fourier differentiation matrix of the specified order with numerical type T on the
domain `x = LinRange{T}(0, 2π, N+1)[1:end-1]`.
"""
function fourier_diff(T::Type{<:Number}, N::Integer; order = 1)
    D = zeros(T, N, N)
    x = LinRange{T}(0, π, N + 1)
    if order == 1
        for i in 2:N
            sgn = (one(T) / 2 - iseven(i))
            D[i, 1] = iseven(N) ? sgn * cot(x[i]) : sgn * csc(x[i])
        end
    elseif order == 2
        D[1, 1] = iseven(N) ? -N^2 * one(T) / 12 - one(T) / 6 :
                  -N^2 * one(T) / 12 + one(T) / 12
        for i in 2:N
            sgn = -(one(T) / 2 - iseven(i))
            D[i, 1] = iseven(N) ? sgn * csc(x[i]) .^ 2 : sgn * cot(x[i]) * csc(x[i])
        end
    else
        error("Not implemented")
    end
    for j in 2:N
        D[1, j] = D[N, j - 1]
        D[2:N, j] .= D[1:(N - 1), j - 1]
    end
    return D
end
fourier_diff(N::Integer; kwargs...) = fourier_diff(Float64, N; kwargs...)

@testitem "Fourier differentiation" begin
    n = 50
    D = NumericalContinuation.fourier_diff(n)
    D2 = NumericalContinuation.fourier_diff(n; order = 2)
    @test_throws ErrorException D3=NumericalContinuation.fourier_diff(n; order = 3)  # for test coverage...
    t = range(0, 2π, n + 1)[1:(end - 1)]
    x = @. sin(2 * cos(t))
    dx_exact = @. -2 * cos(2 * cos(t)) * sin(t)
    dx2_exact = @. -2 * (cos(t) * cos(2 * cos(t)) + 2 * sin(t)^2 * sin(2 * cos(t)))
    @test D * x ≈ dx_exact
    @test D2 * x ≈ dx2_exact
end

struct FourierCollocation{F, T}
    f::F
    u0::Any
    vars::Int
    eqns::Int
    ndim::Int
    nmesh::Int
    Dt::Matrix{T}
    function FourierCollocation(f, u, tspan, p = (), eqns = missing)
        length(tspan) == 2 || throw(ArgumentError("tspan must contain two values"))
        # If needed, call the function to see how many outputs it returns
        ndim = ismissing(eqns) ? length(f(u[:, begin], p, tspan[begin])) : eqns
        _f = IIPWrapper(f, 3)  # Ensure that the function has an in-place form
        T = eltype(u)
        nmesh = size(u, 2)
        Dt = fourier_diff(T, nmesh) .* -2π  # Scale to [0, 1] and transpose (negative = transpose)
        _vars = nmesh * ndim + length(p) + 2
        _eqns = nmesh * ndim
        return new{typeof(_f), T}(_f,
                                  (u = vec(u), p = p, tspan = [tspan[begin], tspan[end]]),
                                  _vars, _eqns, ndim, nmesh, Dt)
    end
end

function (fourier::FourierCollocation{F, TT})(res, uu, data; kwargs...) where {F, TT}
    u = uu.zero
    # TODO: work out if there are any allocations left in here
    # Calculate the right-hand side
    T = u.tspan[end] - u.tspan[begin]
    t = range(u.tspan[begin], u.tspan[end], fourier.nmesh + 1)  # plus one because the final time value is not stored in the state vector (Fourier assumes periodicity)
    idx = 1:(fourier.ndim)
    for i in 1:(fourier.nmesh)
        fourier.f(view(res, idx), view(u.u, idx), u.p, t[i])
        @views res[idx] .*= T
        idx = idx .+ fourier.ndim
    end
    # Compute the difference of the time derivatives using a Fourier differentiation matrix
    u_mat = reshape(u.u, (fourier.ndim, fourier.nmesh))
    res_mat = reshape(res, (fourier.ndim, fourier.nmesh))
    mul!(res_mat, u_mat, fourier.Dt, -one(TT), one(TT))
    return res
end

# Integral phase condition: Int(u_k(t) * u'_{k-1}(t), t=0..1); some codes (e.g., COCO) take
# u_{k-1} to be the initial solution point and don't bother updating. We do the same.
struct PhaseCondition{U}
    du::U  # TODO: think about how to put this into chart data to avoid problems with recomputation
end

# Fourier integration matrix is simply the trapezium rule
(pc::PhaseCondition)(u, p) = dot(pc.du, u)

"""
    fourier_collocation(f, u, tspan, [p]; [t0], [eqns])

Implement a Fourier-based collocation scheme for discretising a periodic orbit. The vector
field `f` is assumed to be of the same form as used in the SciML/DiffEq ecosystem, i.e.,
`f(u, p, t)`, where `u` is the state vector, `p` are the (continuation) parameters, and `t`
is time.

The state vector `u` is assumed to be a matrix; the columns are the state vectors at
different times. The times are assumed to equispaced between `tspan[1]` and `tspan[2]`. The
final state vector at `tspan[2]` is omitted (due to periodicity it is the same as at
`tspan[1]`).

Note that all parameters in `p` are added as (initially inactive) continuation parameters.
If your problem has many parameters that are not likely to be used for continuation, it can
be better to separate them out using a callable `struct` (or a closure) to encapsulate that
data.

By default, the start time `t0` is fixed as zero. This can be changed to be any numerical
value consistent with the problem, or removed entirely by setting `t0 = nothing`.

The number of (first-order) differential equations is specified by `eqns`; by default this
is assumed equal to `length(u[:, begin])` but can be overridden if needed.

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
n = 20  # number of collocation points
t = range(0.0, 2π, length=n+1)[1:end-1]  # omit the last point
u = Matrix(sol(t))
prob = fourier_collocation(hopf!, u, tspan, p)

```
"""
function fourier_collocation(f, u, tspan, p = (); t0 = 0, phase = true, eqns = nothing)
    _eqns = eqns === nothing ? size(u, 1) : eqns
    fcprob = FourierCollocation(f, u, tspan, p, _eqns)
    prob = ContinuationProblem(fcprob)
    # Add continuation parameters (all inactive to start with)
    add_parameters!(prob, :p, keys(p))
    # Fix start time
    if t0 !== nothing
        add_parameter!(prob, :t0, @optic(_.zero.tspan[begin]); value = t0)
    end
    # Add phase condition
    if phase
        # TODO: should this be normalised?
        du = vec((u * fcprob.Dt) .* (1 / fcprob.nmesh))  # time derivative of initial solution scaled by the number of mesh points
        add_monitor_function!(prob, :phase,
                              monitor_function(PhaseCondition(du); active = false,
                                               value = 0))
    end
    return prob
end

@testitem "Fourier collocation" begin
    using LinearAlgebra: norm

    function hopf!(res, u, p, t)
        ss = u[1]^2 + u[2]^2
        res[1] = p[1] * u[1] - u[2] + p[2] * u[1] * ss
        res[2] = u[1] + p[1] * u[2] + p[2] * u[2] * ss
        return res
    end

    # Define initial solution
    p0 = [1.0, -1.0]
    t = range(0, 2π, length = 21)[1:(end - 1)]
    u0 = [sqrt(p0[1]) .* sin.(t) -sqrt(p0[1]) .* cos.(t)]'
    prob = NumericalContinuation.fourier_collocation(hopf!, u0, (0, 2π), p0)

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
    @test length(res) == length(u)

    # Optimised code
    func = NumericalContinuation.ContinuationFunction(prob)

    # Evaluate
    NumericalContinuation.eval_zero_function!(res, func, u, data, active, monitor)
    @test norm(res) < 1e-12
end
