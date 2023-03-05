"""
    fourier_diff([T=Float64,] N; order=1)

Create a Fourier differentiation matrix of the specified order with numerical type T on the
domain `x = LinRange{T}(0, 2π, N+1)[1:end-1]`.
"""
function fourier_diff(T::Type{<:Number}, N::Integer; order = 1)
    D = zeros(T, N, N)
    n1 = (N - 1) ÷ 2
    n2 = N ÷ 2
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
        _f = IIPWrapper(f)  # Ensure that the function has an in-place form
        T = eltype(u)
        nmesh = size(u, 2)
        Dt = fourier_diff(T, nmesh) .* -2π  # Scale to [0, 1] and transpose
        _vars = nmesh * ndim + length(p) + 2
        _eqns = nmesh * ndim
        return new{typeof(_f), T}(_f,
                                  (u = vec(u), p = p, tspan = (tspan[begin], tspan[end])),
                                  _vars, _eqns, ndim, nmesh, Dt, zeros(T, nmesh))
    end
end

function (fourier::FourierCollocation)(res, u, data; parent)
    u_mat = reshape(u.u, (fourier.ndim, fourier.nmesh))
    res_mat = reshape(res, (fourier.ndim, fourier.nmesh))
    # Calculate the right-hand side
    T = u.tspan[end] - u.tspan[begin]
    t = range(u.tspan[begin], u.tspan[end], fourier.nmesh + 1)  # plus one because the final time value is not stored in the state vector (Fourier assumes periodicity)
    for i in Base.OneTo(fourier.nmesh)
        fourier.f(@view(res_mat[:, i]), @view(u_mat[:, i]), u.p, t[i])
        res_mat[:, i] .*= T
    end
    # Compute the difference of the time derivatives using a Fourier differentiation matrix
    mul!(res_mat, u_mat, fourier.Dt, -1, 1)
    return res
end

function fourier_collocation(f, u, tspan, p = (), eqns = missing)
    prob = ContinuationProblem(FourierCollocation(f, u, tspan, p, eqns))
    return add_parameters!(prob, :p, keys(p))
end
