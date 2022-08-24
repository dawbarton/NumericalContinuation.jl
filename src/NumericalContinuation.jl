module NumericalContinuation

using DocStringExtensions

abstract type AbstractContinuationProblem end

struct ContinuationProblem <: AbstractContinuationProblem
    zero_function::Any
    monitor_function::Vector{Any}
    monitor_function_names::Vector{Symbol}
    sub_problem::Vector{Any}
    sub_problem_names::Vector{Symbol}
end

function ContinuationProblem(zero_function)
    ContinuationProblem(zero_function, [], Symbol[], [], Symbol[])
end

function close_problem(prob::ContinuationProblem)
    sub_problem = []
    for i in eachindex(prob.sub_problem)
        push!(sub_problem, close_problem(prob.sub_problem[i]))
    end
    return ClosedProblem(prob.zero_function, prob.monitor_function,
                         prob.monitor_function_names, sub_problem, prob.sub_problem_names)
end

struct ClosedProblem{Z, M, MNAME, P, PNAME} <: AbstractContinuationProblem
    zero_function::Z
    monitor_function::M
    sub_problem::P
end

function ClosedProblem(zero_function, monitor_function, monitor_function_names, sub_problem,
                       sub_problem_names)
    ClosedProblem{typeof(zero_function), Tuple{(typeof.(monitor_function))...},
                  Tuple{monitor_function_names...}, Tuple{(typeof.(sub_problem))...},
                  Tuple{sub_problem_names...}}(zero_function, (monitor_function...,),
                                               (sub_problem...,))
end

# TODO: disallow "zero" and "monitor" as names of subproblems and monitor functions

function zero_function(prob::ContinuationProblem, u, data)
    throw(ArgumentError("Must use close_problem first"))
end

@generated function zero_function(prob::ClosedProblem, u, data)
    return :(ComponentVector($(gen_zero_function(prob, :prob, :u, :data))))
end

"""
    $SIGNATURES

INTERNAL ONLY

Generate an expression tree to call `zero_function` on a `ClosedProblem` and its
subproblems, returning a tuple. This function should be called from an `@generated`
function. `prob`, `u`, and `data` are the names of the corresponding parameters in the
generated function given as a `Symbol` or `Expr`.
"""
function gen_zero_function(::Type{ClosedProblem{Z, M, MNAME, P, PNAME}},
                           prob, u, data) where {Z, M, MNAME, P, PNAME}
    tpl = :(())
    for i in Base.OneTo(length(P.parameters))
        name = PNAME.parameters[i]
        prob_tpl = gen_zero_function(P.parameters[i], :($prob.sub_problem[i]),
                                     :($u.$name),
                                     :($data.$name))
        if prob_tpl != :(())
            push!(tpl.args, :($name = $prob_tpl))
        end
    end
    if Z !== Nothing
        push!(tpl.args,
              :(zero = zero_function($prob.zero_function, $u.zero, $data.zero;
                                     parent = ($u, $data))))
    end
    return tpl
end

function monitor_function(prob::ContinuationProblem, u, data)
    throw(ArgumentError("Must use close_problem first"))
end

@generated function monitor_function(prob::ClosedProblem, u, data)
    return :(ComponentVector($(gen_monitor_function(prob, :prob, :u, :data))))
end

"""
    $SIGNATURES

INTERNAL ONLY

Generate an expression tree to call `monitor_function` on a `ClosedProblem` and its
subproblems, returning a tuple. This function should be called from an `@generated`
function. `prob`, `u`, and `data` are the names of the corresponding parameters in the
generated function given as a `Symbol` or `Expr`.
"""
function gen_monitor_function(::Type{ClosedProblem{Z, M, MNAME, P, PNAME}},
                              prob, u, data) where {Z, M, MNAME, P, PNAME}
    tpl = :(())
    for i in Base.OneTo(length(P.parameters))
        name = PNAME.parameters[i]
        prob_tpl = gen_monitor_function(P.parameters[i], :($prob.sub_problem[i]),
                                        :($u.$name),
                                        :($data.$name))
        if prob_tpl != :(())
            push!(tpl.args, :($name = $prob_tpl))
        end
    end
    for i in Base.OneTo(length(M.parameters))
        name = MNAME.parameters[i]
        push!(tpl.args,
              :($name = monitor_function($prob.monitor_function[i], $u.zero, $data.$name;
                                         parent = ($u, $data))))
    end
    return tpl
end

function monitor_function_names(prob::ContinuationProblem)
    throw(ArgumentError("Must use close_problem first"))
end

function monitor_function_names(prob::ClosedProblem)
    return _monitor_function_names!(Symbol[], typeof(prob), Symbol())
end

function _monitor_function_names!(names::Vector{Symbol},
                                  ::Type{ClosedProblem{Z, M, MNAME, P, PNAME}},
                                  prefix) where {Z, M, MNAME, P, PNAME}
    # Must match the ordering of gen_monitor_function
    for i in Base.OneTo(length(P.parameters))
        name = PNAME.parameters[i]
        _monitor_function_names!(names, P.parameters[i], Symbol(prefix, name, :.))
    end
    for i in Base.OneTo(length(M.parameters))
        name = MNAME.parameters[i]
        push!(names, Symbol(prefix, name))
    end
    return names
end

function test()
    a = ClosedProblem{typeof(sin), Nothing, Nothing, Tuple{}, Tuple{}}(sin, nothing,
                                                                       ())
    b = ClosedProblem{typeof(cos), Nothing, Nothing, Tuple{typeof(a)}, Tuple{:sinprob
                                                                             }}(cos,
                                                                                nothing,
                                                                                (a,))
end

#     res = :(ComponentVector())
#     args = res.args
#     for i in Base.OneTo(n_subproblem)
#         subprob = Symbol(:sub, i)
#         push!(args,
#               Expr(:kw,
#                    subprob,
#                    :(zero_function(prob.sub_problem[$i], u.$subprob, data.$subprob))))
#     end
#     if include_zero
#         push!(args,
#               Expr(:kw, :zero,
#                    :(prob.zero_function(u.zero, data.zero; parent = (u, data)))))
#     end
#     return res
# end

# Handle

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
