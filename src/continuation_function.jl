"""
    count_members(itr)

Count the number of times each element of an iterable occurs. Returns a `Dict` mapping
elements to their count.
"""
function count_members(itr)
    # Taken from https://stackoverflow.com/a/39167771
    d = Dict{eltype(itr), Int}()
    for val in itr
        if isa(val, Number) && isnan(val)
            continue
        end
        d[val] = get(d, val, 0) + 1
    end
    return d
end

"""
    _check_names(names)

INTERNAL ONLY

Check the list of (monitor function and subproblem) names for duplicates and reserved names.
"""
function _check_names(names)
    unique_names = count_members(names)
    dups = findall((>)(1), unique_names)
    if !isempty(dups)
        throw(ArgumentError("Duplicate name within problem: $dups"))
    end
    if any(in.(names, Ref(RESERVED_NAMES)))
        throw(ArgumentError("Reserved names ($RESERVED_NAMES) cannot be used"))
    end
    return
end

"""
    ContinuationFunction

This is a closed (i.e., the structure cannot be modified any further) version of the
`ContinuationProblem` that is used for computational efficiency during continuation.
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
            throw(ArgumentError("Each subproblem of a ContinuationFunction must also be a ContinuationFunction"))
        end
        _check_names([monitor_function_name; sub_problem_name])  # check for duplicates and reserved names
        return new{typeof(zero_function!), Tuple{(typeof.(monitor_function))...},
                   Tuple{monitor_function_name...}, Tuple{(typeof.(sub_problem))...},
                   Tuple{sub_problem_name...}}(zero_function!, (monitor_function...,),
                                               (sub_problem...,))
    end
end

"""
    ContinuationFunction(prob)

Close the `ContinuationProblem` so that structure can no longer be modified (e.g., no new
zero functions / monitor functions can be added). The resulting `ContinuationFunction` can
be used efficiently with functions such as `zero_function` and `monitor_function`.
"""
function ContinuationFunction(prob::ContinuationProblem)
    sub_problem = []
    for i in eachindex(prob.sub_problem)
        push!(sub_problem, ContinuationFunction(prob.sub_problem[i]))
    end
    return ContinuationFunction(prob.zero_function!, prob.monitor_function,
                                prob.monitor_function_name, sub_problem,
                                prob.sub_problem_name)
end

(func::ContinuationFunction)(res, u, (data, chart, active, monitor)) = eval_function!(res, func, u, data, chart, active, monitor)

"""
    FuncPar(name; unwrap=false, append=false)

This struct represents a parameter to a function and is used in the automatic code
generation functions for `ContinuationFunction`s.

- `name` is a symbol or expression passed to the function.
- `unwrap` indicates whether `name` should be unwrapped (via `getproperty`) with the name of
  the subproblem. (E.g., `u` becomes `u.sub_problem_name`.)
- `append` indicates whether `name` should be appended with `zero` or the name of the
  monitor function. (E.g., `u` becomes `u.zero` for zero functions or
  `u.monitor_function_name` for monitor functions.)
"""
struct FuncPar
    name::Any
    unwrap::Bool
    append::Bool
end
FuncPar(name; unwrap=false, append=false) = FuncPar(name, unwrap, append)

"""
    FuncPar(fp::FuncPar, name)

Extend the name of an existing `FuncPar`, respecting the value of `unwrap`.
"""
function FuncPar(fp::FuncPar, name)
    if fp.unwrap
        return FuncPar(:($(fp.name).$name), fp.unwrap, fp.append)
    else
        return fp
    end
end

"""
    get_par(fp::FuncPar, name)

Get the name of a `FuncPar` with the given name appended, respecting the value of `append`.
"""
function get_par(fp::FuncPar, name)
    if fp.append
        return :($(fp.name).$name)
    else
        return fp.name
    end
end

function _generate_apply_function(func::Type{CF}, applyto, args) where {CF <: ContinuationFunction}
    expr = []
    _generate_apply_function!(expr, func, applyto, args)
    return expr
end

function _generate_apply_function!(expr, func::Type{CF}, applyto, args) where {CF <: ContinuationFunction}
    # ContinuationFunction{Z, M, MNAME, P, PNAME}
    (Z, M, MNAME, P, PNAME) = func.parameters
    for i in Base.OneTo(length(P.parameters))
        name = PNAME.parameters[i]
        new_args = [FuncPar(arg, name) for arg in args]
        _generate_apply_function!(expr, P.parameters[i], applyto, new_args)
    end
    if (Z !== Nothing) && ((applyto == :zero) || (applyto == :all))
        apply = Expr(:call)
        for arg in args
            push!(apply.args, get_par(arg, :zero))
        end
        push!(expr, apply)
    end
    if (applyto == :monitor) || (applyto == :all)
        for i in Base.OneTo(length(M.parameters))
            name = MNAME.parameters[i]
            apply = Expr(:call)
            for arg in args
                push!(apply.args, get_par(arg, name))
            end
            push!(expr, apply)
        end
    end
    return expr
end

@generated function eval_function!(res, func::ContinuationFunction, u, data, chart, active,
                                   monitor)
    return _eval_function!(res, func, u, data, chart, active, monitor)
end

function _eval_function!(res, func, u, data, chart, active, monitor)
    result = quote
        j = 0
    end
    _gen_zero_function!(result.args, func, :res, :func, :u, :data, :chart)
    allmonitor_function = _gen_monitor_function!([], func, :func, :u, :data, :chart)
    for (i, monitor_function) in enumerate(allmonitor_function)
        push!(result.args,
              :(monitor_val = active[$i] ? u.monitor[j += 1] : monitor[$i])) # could replace with ifelse and get (with default value);
        # :(monitor_val = ifelse(active[$i], get(u.monitor, j += active[$i], zero(eltype(monitor))), monitor[$i])))  # almost works (doesn't work if u.monitor is missing, e.g., NamedTuple); TODO: benchmark any differences
        push!(result.args,
              :(res.monitor[$i] = $monitor_function - monitor_val))
    end
    push!(result.args, :(return res))
    return result
end

function _eval_function!(res, func::ContinuationFunction, u, data, chart, active, monitor)
    return _eval_function!(typeof(res), typeof(func), typeof(u), typeof(data),
                           typeof(chart), typeof(active), typeof(monitor))
end

@generated function eval_monitor_function!(res, func::ContinuationFunction, u, data, chart)
    result = quote end
    allmonitor_function = _gen_monitor_function!([], func, :func, :u, :data, :chart)
    for (i, monitor_function) in enumerate(allmonitor_function)
        push!(result.args, :(res[$i] = $monitor_function))
    end
    push!(result.args, :(return res))
    return result
end

"""
    _gen_zero_function!(result, ::ContinuationFunction, res, func, u, data, chart)

INTERNAL ONLY

Generate an expression tree to call each zero function within a `ContinuationFunction` and
its subproblems, acting in place. This function should be called from an `@generated`
function. `res`, `prob`, `u`, `data`, `chart` are the names of the corresponding parameters
in the generated function given as a `Symbol` or `Expr`.
"""
function _gen_zero_function!(result, ::Type{ContinuationFunction{Z, M, MNAME, P, PNAME}},
                             res, func, u, data, chart) where {Z, M, MNAME, P, PNAME}
    for i in Base.OneTo(length(P.parameters))
        name = PNAME.parameters[i]
        _gen_zero_function!(result, P.parameters[i], :($res.$name),
                            :($func.sub_problem[$i]), :($u.$name), :($data.$name), chart)
    end
    if Z !== Nothing
        push!(result,
              :($func.zero_function!($res.zero, $u.zero, $data.zero; parent = ($u, $data),
                                     chart = $chart)))
    end
    return result
end

"""
    _gen_monitor_function!(result, ::ContinuationFunction, res, prob, u, data, chart)

INTERNAL ONLY

Generate an expression tree to call each monitor function within a `ContinuationFunction`
and its subproblems, acting in place. This function should be called from an `@generated`
function. `res`, `func`, `u`, `data`, and `chart` are the names of the corresponding
parameters in the generated function given as a `Symbol` or `Expr`.
"""
function _gen_monitor_function!(result, ::Type{ContinuationFunction{Z, M, MNAME, P, PNAME}},
                                func, u, data, chart) where {Z, M, MNAME, P, PNAME}
    for i in Base.OneTo(length(P.parameters))
        name = PNAME.parameters[i]
        _gen_monitor_function!(result, P.parameters[i], :($func.sub_problem[$i]),
                               :($u.$name), :($data.$name), chart)
    end
    for i in Base.OneTo(length(M.parameters))
        name = MNAME.parameters[i]
        push!(result,
              :($func.monitor_function[$i]($u.zero, $data.$name; parent = ($u, $data),
                                           chart = $chart)))
    end
    return result
end
