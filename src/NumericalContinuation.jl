module NumericalContinuation

# TODO: the _.zero I added to add_parameter! is going to conflict with the different calling
# conventions of monitor_function (i.e., f(u) versus f(u, p)); maybe shift the parameters to
# add_parameter! as a closer-to-the-user function

using TestItems: @testitem
using LinearAlgebra
using ComponentArrays: ComponentVector, getaxes
using Accessors: @optic, opcompose, PropertyLens, IndexLens
using ForwardDiff: jacobian, jacobian!, JacobianConfig
using SciMLBase: SciMLBase, solve
using NonlinearSolve: NonlinearSolve, NonlinearProblem, init, reinit!, solve!,
                      NewtonRaphson, TrustRegion, ReturnCode
using PreallocationTools: DiffCache, get_tmp
using FFTW: ifft!
using QuadGK: quadgk, gauss

const RESERVED_NAMES = Set([:zero, :monitor])

export solve, PseudoArclength

include("continuation_problem.jl")
include("continuation_function.jl")
include("algebraic_problems.jl")
include("polynomials.jl")
include("polynomial_collocation.jl")
include("fourier_collocation.jl")
include("covering.jl")

export test_problem0

function test_problem0(T = Float64)
    return algebraic_problem((u, p) -> [u[1]^3 - u[1] - p[1]], T[1], T[0])
end

export test_problem1

function test_problem1(T = Float64)
    _u0 = T[1.0, 1.0, 0.5]
    prob = algebraic_problem((u, p) -> [u[1]^2 + (u[2] - 1)^2 - 1], _u0)
    add_monitor_function!(prob, :ψ₁, monitor_function((u, p) -> sqrt(u[1]^2 + u[2]^2)))
    add_monitor_function!(prob, :ψ₂, monitor_function((u, p) -> u[1]))
    add_monitor_function!(prob, :ψ₃, monitor_function((u, p) -> u[3] - u[1] + 0.5))
    prob2 = algebraic_problem((u, p) -> [u[1]^2 + (u[2] - 1)^2 - 1], _u0)
    add_monitor_function!(prob2, :ψ₁, monitor_function((u, p) -> sqrt(u[1]^2 + u[2]^2)))
    add_monitor_function!(prob2, :ψ₂, monitor_function((u, p) -> u[1]))
    add_monitor_function!(prob2, :ψ₃, monitor_function((u, p) -> u[3] - u[1] + 0.5))
    add_sub_problem!(prob, :subprob, prob2)
    return prob
end

export test_problem2

function test_problem2(T = Float64)
    function hopf!(res, u, p, t)
        ss = u[1]^2 + u[2]^2
        res[1] = p[1] * u[1] - u[2] + p[2] * u[1] * ss
        res[2] = u[1] + p[1] * u[2] + p[2] * u[2] * ss
        return res
    end

    # Define initial solution
    p0 = [1.0, -1.0]
    t = range(0, 2π, length = 21)[1:(end - 1)]
    u0 = 0.9 .* [sqrt(p0[1]) .* cos.(t) sqrt(p0[1]) .* sin.(t)]'
    prob = NumericalContinuation.fourier_collocation(hopf!, u0, (0, 2π), p0; phase = false)
    # Don't use the integral phase condition; fix u[2] = 0 instead
    add_parameter!(prob, :phase, @optic _.u[2]; value = 0)
    return prob
end

export test_problem3

function duffing!(du, u, p, t)
    du[1] = u[2]
    du[2] = p.Γ * sin(p.ω * t + p.ϕ) - p.ξ * u[2] - p.k * u[1] - p.k₃ * u[1]^3
end

function test_problem3(; nmesh = 20)
    p0 = (Γ = 1.0, ω = 0.1, ϕ = π / 2, ξ = 0.05, k = 1.0, k₃ = 0.1)
    u0 = t -> [cos(t), -p0.ω*sin(t)]
    prob = NumericalContinuation.limit_cycle(duffing!, u0, (0, 2π / p0.ω), p0;
                                                     phase = false, nmesh)
    add_parameter!(prob, :period, u -> u.coll.tspan[2] - 2π / u.coll.p.ω; value = 0)
    add_parameter!(prob, :phase, @optic _.coll.u[2]; value = 0)
    return prob
end

end
