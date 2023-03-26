module NumericalContinuation

using TestItems: @testitem
using LinearAlgebra
using ComponentArrays: ComponentVector, getaxes, label2index
using Accessors: @optic, opcompose, PropertyLens, IndexLens
using ForwardDiff
using NonlinearSolve: NonlinearProblem, solve, TrustRegion, ReturnCode

const RESERVED_NAMES = Set([:zero, :monitor])

include("continuation_problem.jl")
include("continuation_function.jl")
include("algebraic_problems.jl")
include("boundary_value_problems.jl")

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

# prob = zero_problem(u -> [u.x^2 + u.y^2 - 1]; u0=ComponentArray(x=1.0, y=0.0))
# prob = zero_problem(u -> [u[1]^2 + u[2]^2 - 1]; u0=[1.0, 0.0])
# prob = periodic_orbit_problem(u -> [u[2], -u[1]]; u0=[1.0, 2.0], T=2.5)

function trace_branch(prob)
    chart = create_initial_chart(prob)
    branch = create_branch(prob, chart)
    ctr = 1
    max_steps = get_option(prob, :max_steps)
    while ctr < max_steps
        pred = predict_from_chart(prob, chart)
        chart = correct_chart(prob, chart)
        push!(branch, chart)
        ctr += 1
    end
end

end
