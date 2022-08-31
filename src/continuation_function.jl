
"""
    $TYPEDEF

This is a closed (i.e., the structure cannot be modified any further) version of the
`ContinuationProblem` that is used for computational efficiency during continuation.

See [`close_problem`](@ref)
"""
struct ContinuationFunction{Z, M, MNAME, P, PNAME}
    zero_function!::Z
    monitor_function::M
    sub_problem::P

    function ContinuationFunction(zero_function!, monitor_function, monitor_function_name,
                                  sub_problem, sub_problem_name)
        if length(monitor_function) != length(monitor_function_name)
            throw(ArgumentError("Each monitor function must have one and only one name"))
        end
        if length(sub_problem) != length(sub_problem_name)
            throw(ArgumentError("Each subproblem must have one and only one name"))
        end
        if !all(isa.(sub_problem, ContinuationFunction))
            throw(ArgumentError("All subproblems of a closed problem must also be closed"))
        end
        _check_names([monitor_function_name; sub_problem_name])  # check for duplicates and reserved names
        return new{typeof(zero_function!), Tuple{(typeof.(monitor_function))...},
                   Tuple{monitor_function_name...}, Tuple{(typeof.(sub_problem))...},
                   Tuple{sub_problem_name...}}(zero_function!, (monitor_function...,),
                                               (sub_problem...,))
    end
end

"""
    $SIGNATURES

Close the `ContinuationProblem` so that structure can no longer be modified (e.g., no new
zero functions / monitor functions can be added). The resulting `ContinuationFunction` can be used
efficiently with functions such as `zero_function` and `monitor_function`.
"""
function ContinuationFunction(prob::ContinuationProblem)
    sub_problem = []
    for i in eachindex(prob.sub_problem)
        push!(sub_problem, close_problem(prob.sub_problem[i]))
    end
    return ContinuationFunction(prob.zero_function!, prob.monitor_function,
                                prob.monitor_function_name, sub_problem,
                                prob.sub_problem_name)
end

@generated function eval_function!(res, prob::ContinuationFunction, u, data, active,
                                   monitor)
    result = quote
        j = 0
    end
    _gen_zero_function!(result.args, prob, :(res.zero), :prob, :u, :data)
    allmonitor_function = _gen_monitor_function!([], prob, :(res.monitor), :prob, :u, :data)
    for (i, monitor_function) in enumerate(allmonitor_function)
        push!(result.args,
              :(monitor_val = $i in active ? u.monitor[j += 1] : monitor[$i])) # could replace with ifelse and get (with default value)
        push!(result.args,
              :(res.monitor[$i] = $monitor_function - monitor_val))
    end
    push!(result.args, :(return res))
    return result
end

@generated function eval_monitor_function!(res, prob::ContinuationFunction, u, data)
    result = quote end
    allmonitor_function = _gen_monitor_function!([], prob, :(res.monitor), :prob, :u, :data)
    for (i, monitor_function) in enumerate(allmonitor_function)
        push!(result.args, :(res[$i] = $monitor_function))
    end
    push!(result.args, :(return res))
    return result
end

"""
    $SIGNATURES

INTERNAL ONLY

Generate an expression tree to call each zero function within a `ContinuationFunction` and its
subproblems, acting in place. This function should be called from an `@generated`
function. `prob`, `u`, and `data` are the names of the corresponding parameters in the
generated function given as a `Symbol` or `Expr`.
"""
function _gen_zero_function!(result, ::Type{ContinuationFunction{Z, M, MNAME, P, PNAME}},
                             res, prob, u, data) where {Z, M, MNAME, P, PNAME}
    for i in Base.OneTo(length(P.parameters))
        name = PNAME.parameters[i]
        _gen_zero_function!(result, P.parameters[i], :($res.$name),
                            :($prob.sub_problem[$i]), :($u.$name), :($data.$name))
    end
    if Z !== Nothing
        push!(result,
              :($prob.zero_function!($res.zero, $u.zero, $data.zero; parent = ($u, $data))))
    end
    return result
end

"""
    $SIGNATURES

INTERNAL ONLY

Generate an expression tree to call each monitor function within a `ContinuationFunction` and its
subproblems, acting in place. This function should be called from an `@generated`
function. `prob`, `u`, and `data` are the names of the corresponding parameters in the
generated function given as a `Symbol` or `Expr`.
"""
function _gen_monitor_function!(result, ::Type{ContinuationFunction{Z, M, MNAME, P, PNAME}},
                                res, prob, u, data) where {Z, M, MNAME, P, PNAME}
    for i in Base.OneTo(length(P.parameters))
        name = PNAME.parameters[i]
        _gen_monitor_function!(result, P.parameters[i], :($res.$name),
                               :($prob.sub_problem[$i]), :($u.$name), :($data.$name))
    end
    for i in Base.OneTo(length(M.parameters))
        name = MNAME.parameters[i]
        push!(result,
              :($prob.monitor_function[$i]($u.zero, $data.$name; parent = ($u, $data))))
    end
    return result
end
