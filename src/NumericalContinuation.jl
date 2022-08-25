module NumericalContinuation

using ComponentArrays
using DocStringExtensions

const RESERVED_NAMES = Set([:zero])

include("continuation_problem.jl")

function test()
    a = ContinuationProblem(sin, monitor=[:a=>identity, :b=>identity])
    return ContinuationProblem(cos, sub = [:sin => a])
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
