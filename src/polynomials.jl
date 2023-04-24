# Barycentric formula for polynomial interpolation on equispaced points, Chebyshev points of
# the second kind, and Legendre points. The formulae used are taken from the paper of Berrut
# and Trefethen, SIAM Review, 2004.
#
# Also see: Chebfun (legpts.m, chebpts.m)

"""
    OrthogonalPolynomial{T}

Represents an orthogonal polynomial with underlying numerical type `T`. Enables fast
interpolation, differentiation, and integration. The domain is generally taken as [-1, +1].
"""
struct OrthogonalPolynomial{T}
    "Nodes of the polynomial"
    nodes::Vector{T}
    "Weights for the barycentric formula"
    bary_weights::Vector{T}
    "Weights for the quadrature formula"
    quad_weights::Vector{T}
    "Indicates whether the domain is closed"
    closed::Bool
end

"""
    nodes(poly::OrthogonalPolynomial)

Returns the nodes of the polynomial.
"""
nodes(poly::OrthogonalPolynomial) = poly.nodes

"""
    barycentric_weights(poly::OrthogonalPolynomial)

Returns the weights for the barycentric interpolation formula.
"""
barycentric_weights(poly::OrthogonalPolynomial) = poly.bary_weights

"""
    quadrature_weights(poly::OrthogonalPolynomial)

Returns the weights for the corresponding quadrature formula.
"""
quadrature_weights(poly::OrthogonalPolynomial) = poly.quad_weights

"""
    is_closed_interval(poly::OrthogonalPolynomial)

Returns `true` if the domain is closed (i.e., the end points of the domain are a node of the
polynomical), `false` otherwise.
"""
is_closed_interval(poly::OrthogonalPolynomial) = poly.closed

function chebyshev_nodes(T, n::Integer)
    return T[-cospi(T(j) / n) for j in 0:n]
end

function chebyshev_barycentric_weights(T, n::Integer)
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

function chebyshev_quadrature_weights(T, n::Integer)
    if n == 0
        return T[2]
    elseif n > 0
        c = Vector{Complex{T}}(undef, n + 1)
        w = Vector{T}(undef, n + 1)
        for i in 1:(n ÷ 2 + 1)
            c[i] = T(2) / T(1 - (2 * (i - 1))^2)
        end
        for i in (n ÷ 2 + 2):n
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

"""
    chebyshev([T], n::Integer)

Return an `OrthogonalPolynomial` representing the Chebyshev polynomial of the second kind
with order `n`.
"""
function chebyshev(T, n::Integer)
    return OrthogonalPolynomial{T}(chebyshev_nodes(T, n),
                                   chebyshev_barycentric_weights(T, n),
                                   chebyshev_quadrature_weights(T, n),
                                   true)
end

chebyshev(n::Integer) = chebyshev(Float64, n)

@testitem "Chebyshev polynomials" begin
    testfcn = x -> sinpi(x) * exp(x)
    dtestfcn = x -> (π * cospi(x) + sinpi(x)) * exp(x)
    for n in (10, 11)  # test both odd and even polynomial orders
        cheb = NumericalContinuation.chebyshev(n)
        x0 = NumericalContinuation.nodes(cheb)
        y0 = testfcn.(x0)
        x = range(-1, stop = 1, length = 101)
        y = NumericalContinuation.interpolate(cheb, y0, x)
        @test y≈testfcn.(x) atol=1e-3
        In = NumericalContinuation.interpolation_matrix(cheb, x)
        @test In * y0 ≈ y
        D = NumericalContinuation.differentiation_matrix(cheb)
        dy0 = D * y0
        @test dy0≈dtestfcn.(x0) atol=1e-3
        int = NumericalContinuation.integrate(cheb, dtestfcn.(x0))
        @test int≈testfcn(1) - testfcn(-1) atol=1e-4
    end
end

function equispaced_nodes(T, n::Integer)
    return collect(LinRange{T}(-1, 1, n + 1))
end

function equispaced_barycentric_weights(T, n::Integer)
    return T[((2 * xor(isodd(n), iseven(j)) - 1) * binomial(n, j)) for j in 0:n]
end

# Newton-Cotes coefficients taken from https://oeis.org/A093735 and https://oeis.org/A093736
const NEWTON_COTES_NUM = ([1, 1], [1, 4, 1], [3, 9, 9, 3], [14, 64, 8, 64, 14],
                          [95, 125, 125, 125, 125, 95], [41, 54, 27, 68, 27, 54, 41],
                          [5257, 25039, 343, 20923, 20923, 343, 25039, 5257],
                          [3956, 23552, -3712, 41984, -3632, 41984, -3712, 23552, 3956])
const NEWTON_COTES_DENOM = ([2, 2], [3, 3, 3], [8, 8, 8, 8], [45, 45, 15, 45, 45],
                            [288, 96, 144, 144, 96, 288], [140, 35, 140, 35, 140, 35, 140],
                            [17280, 17280, 640, 17280, 17280, 640, 17280, 17280],
                            [14175, 14175, 14175, 14175, 2835, 14175, 14175, 14175, 14175])

function equispaced_quadrature_weights(T, n::Integer)
    h = T(2) / n
    return h .* T.(NEWTON_COTES_NUM[n]) ./ T.(NEWTON_COTES_DENOM[n])
end

"""
    equispaced([T], n::Integer)

Return an `OrthogonalPolynomial` representing the a Lagrange polynomial with equispaced
nodes of order `n`.
"""
function equispaced(T, n::Integer)
    return OrthogonalPolynomial{T}(equispaced_nodes(T, n),
                                   equispaced_barycentric_weights(T, n),
                                   equispaced_quadrature_weights(T, n),
                                   true)
end

equispaced(n::Integer) = equispaced(Float64, n)

const lagrange = equispaced

@testitem "Lagrange (equispaced) polynomials" begin
    testfcn = x -> sinpi(x) * exp(x)
    dtestfcn = x -> (π * cospi(x) + sinpi(x)) * exp(x)
    for n in (7, 8)  # test both odd and even polynomial orders
        equi = NumericalContinuation.equispaced(n)
        x0 = NumericalContinuation.nodes(equi)
        y0 = testfcn.(x0)
        x = range(-1, stop = 1, length = 101)
        y = NumericalContinuation.interpolate(equi, y0, x)
        @test y≈testfcn.(x) atol=1e-1
        In = NumericalContinuation.interpolation_matrix(equi, x)
        @test In * y0 ≈ y
        D = NumericalContinuation.differentiation_matrix(equi)
        dy0 = D * y0
        @test dy0≈dtestfcn.(x0) atol=5e-1
        int = NumericalContinuation.integrate(equi, dtestfcn.(x0))
        @test int≈testfcn(1) - testfcn(-1) atol=1e-2
    end
end

"""
    legendre([T], n::Integer)

Return an `OrthogonalPolynomial` representing the Legendre polynomials of order `n`.
"""
function legendre(T, n::Integer)
    (x, w) = gauss(T, n + 1)
    l = T[(2 * isodd(i) - 1) * sqrt((1 - x[i]^2) * w[i]) for i in eachindex(x)]
    return OrthogonalPolynomial{T}(x, l, w, false)
end

legendre(n::Integer) = legendre(Float64, n)

@testitem "Legendre polynomials" begin
    testfcn = x -> sinpi(x) * exp(x)
    dtestfcn = x -> (π * cospi(x) + sinpi(x)) * exp(x)
    for n in (10, 11)  # test both odd and even polynomial orders
        cheb = NumericalContinuation.legendre(n)
        x0 = NumericalContinuation.nodes(cheb)
        y0 = testfcn.(x0)
        x = range(-1, stop = 1, length = 101)
        y = NumericalContinuation.interpolate(cheb, y0, x)
        @test y≈testfcn.(x) atol=1e-3
        In = NumericalContinuation.interpolation_matrix(cheb, x)
        @test In * y0 ≈ y
        D = NumericalContinuation.differentiation_matrix(cheb)
        dy0 = D * y0
        @test dy0≈dtestfcn.(x0) atol=1e-2
        int = NumericalContinuation.integrate(cheb, dtestfcn.(x0))
        @test int≈testfcn(1) - testfcn(-1) atol=1e-4
    end
end

"""
    interpolation_matrix(poly::OrthogonalPolynomial, x)

Return the interpolation matrix from the nodes of `poly` to the point(s) `x`.
For example :

    P = chebyshev(5)
    x = range(-1, stop=1, length=10)
    M = interpolation_matrix(P, x)

Now `y(x) ≈ M*y₀` given that `y(nodes(poly)) = y₀.
"""
function interpolation_matrix(poly::OrthogonalPolynomial{T}, x) where {T}
    w = barycentric_weights(poly)
    x₀ = nodes(poly)
    N1 = length(x₀) # N + 1
    M = Matrix{T}(undef, length(x), N1)
    # Eq. (4.2)
    for j in eachindex(x)
        xx = convert(T, x[j])
        Msum = zero(T)
        exact = 0
        for i in Base.OneTo(N1)
            exact = ifelse(xx == x₀[i], i, exact)
            M[j, i] = w[i] / (xx - x₀[i])
            Msum += M[j, i]
        end
        if Msum == 0
            for i in Base.OneTo(N1)
                M[j, i] = zero(T)
            end
        elseif exact > 0
            for i in Base.OneTo(N1)
                M[j, i] = zero(T)
            end
            M[j, exact] = one(T)
        else
            for i in Base.OneTo(N1)
                M[j, i] /= Msum
            end
        end
    end
    return M
end

"""
    interpolate(poly::OrthogonalPolynomial, y₀, [x])

Return the value of `y(x)` given that `y(nodes(poly)) = y₀`. If the value of `x`
is not provided, return a function `y(x)` that evaluates the interpolant at any
`x`.
"""
function interpolate(poly::OrthogonalPolynomial{T}, y₀, x::Number) where {T}
    w = barycentric_weights(poly)
    x₀ = nodes(poly)
    N1 = length(x₀) # N + 1
    xx = convert(T, x)
    # Eq. (4.2)
    numer = zero(T)
    denom = zero(T)
    exact = 0
    for j in Base.OneTo(N1)
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

function interpolate(poly::OrthogonalPolynomial{T}, y₀, x) where {T}
    return [interpolate(poly, y₀, xᵢ) for xᵢ in x]
end

function interpolate(poly::OrthogonalPolynomial{T}, y₀) where {T}
    return x -> interpolate(poly, y₀, x)
end

"""
    differentiation_matrix(poly::OrthogonalPolynomial)

Return the differentiation matrix at the nodes of the polynomial specified.

    P = chebyshev(5)
    D = differentiation_matrix(P)

Now dy/dx ≈ `D*y` at the nodes of the polynomial.
"""
function differentiation_matrix(poly::OrthogonalPolynomial{T}) where {T}
    # Eqs. (9.4) and (9.5)
    w = barycentric_weights(poly)
    x = nodes(poly)
    N1 = length(x) # N + 1
    D = Matrix{T}(undef, N1, N1)
    for i in Base.OneTo(N1)
        Dsum = zero(T)
        for j in Base.OneTo(i - 1)
            temp = (w[j] / w[i]) / (x[i] - x[j])
            D[i, j] = temp
            Dsum += temp
        end
        for j in (i + 1):(N1)
            temp = (w[j] / w[i]) / (x[i] - x[j])
            D[i, j] = temp
            Dsum += temp
        end
        D[i, i] = -Dsum
    end
    return D
end

"""
    integrate(poly::OrthogonalPolynomial, y₀)

Return the integral of `y₀` over the domain of the polynomial specified (typically [-1, +1];
can be determined using `nodes(poly)`).
"""
integrate(poly::OrthogonalPolynomial, y₀) = sum(dot(quadrature_weights(poly), y₀))
