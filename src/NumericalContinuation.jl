module NumericalContinuation

using TestItems: @testitem
using LinearAlgebra
using ComponentArrays: ComponentVector, getaxes
using Accessors: @optic, opcompose, PropertyLens, IndexLens
using ForwardDiff: jacobian, jacobian!, JacobianConfig
using SciMLBase: SciMLBase, solve
using NonlinearSolve: NonlinearSolve, NonlinearProblem, init, reinit!, solve!, NewtonRaphson, TrustRegion, ReturnCode
using PreallocationTools: DiffCache, get_tmp

const RESERVED_NAMES = Set([:zero, :monitor])

export solve, PseudoArclength

include("continuation_problem.jl")
include("continuation_function.jl")
include("algebraic_problems.jl")
include("boundary_value_problems.jl")
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
    u0 = 0.9 .* [sqrt(p0[1]) .* sin.(t) -sqrt(p0[1]) .* cos.(t)]'
    return NumericalContinuation.fourier_collocation(hopf!, u0, (0, 2π), p0)
end

end
