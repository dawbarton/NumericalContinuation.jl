# Barycentric formula for polynomial interpolation on equispaced points, Chebyshev points of
# the second kind, and Legendre points. The formulae used are taken from the paper of Berrut
# and Trefethen, SIAM Review, 2004.
#
# Also see: Chebfun (legpts.m, chebpts.m)

struct OrthogonalPolynomial{T}
    "Nodes of the polynomial"
    nodes::Vector{T}
    "Weights for the barycentric formula"
    bary_weights::Vector{T}
    "Weights for the quadrature formula"
    quad_weights::Vector{T}
end

function chebyshev_nodes(T, n)
    return T[-cospi(T(j) / N) for j in 0:n]
end

function chebyshev_barycentric_weights(T, n)
    w = Vector{T}(undef, n + 1)
    w[1] = -0.5
    sgn = +1
    for i in 2:n
        w[i] = sgn
        sgn = -sgn
    end
    w[n + 1] = ifelse(sgn < 0, -0.5, 0.5)
    return w
end

function chebyshev_quadrature_weights(T, n)
    if n == 0
        return T[2]
    elseif n > 0
        c = Vector{Complex{T}}(undef, n + 1)
        w = Vector{T}(undef, n + 1)
        for i in 1:n ÷ 2 + 1
            c[i] = T(2) / T(1 - (2 * (i - 1))^2)
        end
        for i in n ÷ 2 + 2:n
            c[i] = c[n - i + 2]
        end
        ifft!(@view(c[1:n]))
        w[1] = real(c[1]) / T(2)
        for i in 2:n
            w[i] = real(c[i])
        end
        w[n + 1] = w[1]
        return w
    else
        error("n must be non-negative")
    end
end

function chebyshev(T, n)

end


function equi_bary_weights(T, N)
    T[((2 * xor(isodd(N), iseven(j)) - 1) * binomial(N, j)) for j in 0:N]
end

# Newton-Cotes coefficients taken from https://oeis.org/A093735 and https://oeis.org/A093736
const NEWTON_COTES_NUM = ((1, 1), (1, 4, 1), (3, 9, 9, 3), (14, 64, 8, 64, 14),
                          (95, 125, 125, 125, 125, 95), (41, 54, 27, 68, 27, 54, 41),
                          (5257, 25039, 343, 20923, 20923, 343, 25039, 5257),
                          (3956, 23552, -3712, 41984, -3632, 41984, -3712, 23552, 3956))
const NEWTON_COTES_DENOM = ((2, 2), (3, 3, 3), (8, 8, 8, 8), (45, 45, 15, 45, 45),
                            (288, 96, 144, 144, 96, 288), (140, 35, 140, 35, 140, 35, 140),
                            (17280, 17280, 640, 17280, 17280, 640, 17280, 17280),
                            (14175, 14175, 14175, 14175, 2835, 14175, 14175, 14175, 14175))

equi_nodes(T, N) = T[(-1 + 2j / N) for j in 0:N]

function nodes(poly::Type{<:Equispaced{N}}, shift::Number, scale::Number) where {N}
    range(shift - scale, stop = shift + scale, length = N + 1)
end

@inline function _node(poly::Type{<:Chebyshev1{N, T}}, j::Integer) where {N, T}
    -cospi(T(2j + 1) / (2N + 2))
end

@inline _node(poly::Type{<:Chebyshev2{N, T}}, j::Integer) where {N, T} = -cospi(T(j) / N)

"""
    interpolation_matrix(poly::AbstractPolynomial, x)

Return the interpolation matrix from the nodes of `poly` to the point(s) `x`.
For example :

    P = Chebyshev2{5}()
    x = range(-1, stop=1, length=10)
    M = interpolation_matrix(P, x)

Now `y(x) ≈ M*y₀` given that `y(nodes(poly)) = y₀.
"""
function interpolation_matrix(poly::AbstractPolynomial{N, T},
                              x::Union{Number, AbstractVector}) where {N, T}
    w = weights(poly)
    x₀ = nodes(poly)
    M = Matrix{T}(undef, length(x), N + 1)
    # Eq. (4.2)
    for j in eachindex(x)
        xx = convert(T, x[j])
        Msum = zero(T)
        exact = 0
        for i in Base.OneTo(N + 1)
            exact = ifelse(xx == x₀[i], i, exact)
            M[j, i] = w[i] / (xx - x₀[i])
            Msum += M[j, i]
        end
        if Msum == 0
            for i in Base.OneTo(N + 1)
                M[j, i] = zero(T)
            end
        elseif exact > 0
            for i in Base.OneTo(N + 1)
                M[j, i] = zero(T)
            end
            M[j, exact] = one(T)
        else
            for i in Base.OneTo(N + 1)
                M[j, i] /= Msum
            end
        end
    end
    return M
end

"""
    interpolate(poly::AbstractPolynomial, y₀, [x])

Return the value of `y(x)` given that `y(nodes(poly)) = y₀`. If the value of `x`
is not provided, return a function `y(x)` that evaluates the interpolant at any
`x`.
"""
function interpolate(poly::AbstractPolynomial{N, T}, y₀::AbstractVector{<:Number},
                     x::Number) where {N, T}
    w = weights(poly)
    x₀ = nodes(poly)
    xx = convert(T, x)
    # Eq. (4.2)
    numer = zero(T)
    denom = zero(T)
    exact = 0
    for j in Base.OneTo(N + 1)
        xdiff = xx - x₀[j]
        temp = w[j] / xdiff
        numer += temp * convert(T, y₀[j])
        denom += temp
        exact = ifelse(xdiff == 0, j, exact)
    end
    if exact > 0
        return convert(T, y₀[exact])
    else
        return numer / denom
    end
end

function interpolate(poly::AbstractPolynomial{N, T},
                     y₀::AbstractVector{<:Number}) where {N, T}
    x -> interpolate(poly, y₀, x)
end
function interpolate(poly::AbstractPolynomial{N, T}, y₀::AbstractVector{<:Number},
                     x::AbstractVector{<:Number}) where {N, T}
    [interpolate(poly, y₀, xᵢ) for xᵢ in x]
end

"""
    differentiation_matrix(poly::AbstractPolynomial)

Return the differentiation matrix at the nodes of the polynomial specified.

    P = Chebyshev2(5)
    D = differentiation_matrix(P)

Now dy/dx ≈ `D*y` at the nodes of the polynomial.
"""
function differentiation_matrix(poly::AbstractPolynomial{N, T}) where {N, T}
    # Eqs. (9.4) and (9.5)
    w = weights(poly)
    x = nodes(poly)
    D = Matrix{T}(undef, N + 1, N + 1)
    for i in Base.OneTo(N + 1)
        Dsum = zero(T)
        for j in Base.OneTo(i - 1)
            temp = (w[j] / w[i]) / (x[i] - x[j])
            D[i, j] = temp
            Dsum += temp
        end
        for j in (i + 1):(N + 1)
            temp = (w[j] / w[i]) / (x[i] - x[j])
            D[i, j] = temp
            Dsum += temp
        end
        D[i, i] = -Dsum
    end
    return D
end
