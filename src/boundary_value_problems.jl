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
        Dt = fourier_diff(T, nmesh) .* -2π  # Scale to [0, 1] and transpose
        _vars = nmesh * ndim + length(p) + 2
        _eqns = nmesh * ndim
        return new{typeof(_f), T}(_f,
                                  (u = vec(u), p = p, tspan = [tspan[begin], tspan[end]]),
                                  _vars, _eqns, ndim, nmesh, Dt)
    end
end

function (fourier::FourierCollocation)(res, u, data; kwargs...)
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
    mul!(res_mat, u_mat, fourier.Dt, -1.0, 1.0)
    return res
end

function fourier_collocation(f, u, tspan, p = (), eqns = missing)
    prob = ContinuationProblem(FourierCollocation(f, u, tspan, p, eqns))
    return add_parameters!(prob, :p, keys(p))
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
    u0 = [sqrt(p0[1]) .* cos.(t) sqrt(p0[1]) .* sin.(t)]'
    prob = NumericalContinuation.fourier_collocation(hopf!, u0, (0, 2π), p0, 2)

    # Problem set up
    (_u, data) = NumericalContinuation.get_initial(prob)
    (_monitor, _active) = NumericalContinuation.get_initial_monitor(prob, _u, data)
    res_layout = NumericalContinuation.get_initial_residual_layout(prob)
    chart = nothing

    # Buffers
    u = ComponentVector{Float64}(_u)
    monitor = ComponentVector{Float64}(_monitor)
    res = ComponentVector{Float64}(res_layout)
    active = ComponentVector{Bool}(_active)

    # Optimised code
    func = NumericalContinuation.ContinuationFunction(prob)

    # Evaluate
    NumericalContinuation.eval_function!(res, func, u, data, chart, active, monitor)
    @test norm(res) < 1e-12
end
