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
    prob::ContinuationProblem
    zero_function!::Z
    monitor_function::M
    sub_problem::P

    function ContinuationFunction(prob, zero_function!, monitor_function,
                                  monitor_function_name, sub_problem, sub_problem_name)
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
                   Tuple{sub_problem_name...}}(prob, zero_function!, (monitor_function...,),
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
    return ContinuationFunction(prob, prob.zero_function!, prob.monitor_function,
                                prob.monitor_function_name, sub_problem,
                                prob.sub_problem_name)
end

"""
    FuncPar(name; unwrap=false, append=false, cont_func=false)

This struct represents a parameter to a function and is used in the automatic code
generation functions for `ContinuationFunction`s.

- `name` is a symbol or expression passed to the function.
- `unwrap` indicates whether `name` should be unwrapped (via `getproperty`) with the name of
  the subproblem. (E.g., `u` becomes `u.sub_problem_name`.)
- `append` indicates whether `name` should be appended with `zero` or the name of the
  monitor function. (E.g., `u` becomes `u.zero` for zero functions or
  `u.monitor_function_name` for monitor functions.)
- `cont_func` indicates that the variable passed is a `ContinuationFunction`.
"""
struct FuncPar
    name::Any
    unwrap::Bool
    append::Bool
    cont_func::Bool
end
function FuncPar(name; unwrap = false, append = false, cont_func = false)
    FuncPar(name, unwrap, append, cont_func)
end

"""
    FuncPar(fp::FuncPar, name, idx)

Extend the name of an existing `FuncPar`, respecting the value of `unwrap`.
"""
function FuncPar(fp::FuncPar, name, idx)
    if fp.cont_func
        return FuncPar(:($(fp.name).sub_problem[$idx]), fp.unwrap, fp.append, fp.cont_func)
    elseif fp.unwrap
        return FuncPar(:($(fp.name).$name), fp.unwrap, fp.append, fp.cont_func)
    else
        return fp
    end
end

"""
    get_par(fp::FuncPar, name, idx)

Get the name of a `FuncPar` with the given name appended, respecting the value of `append`.
"""
function get_par(fp::FuncPar, name, idx)
    if fp.cont_func
        if idx == 0
            return :($(fp.name).zero_function!)
        else
            return :($(fp.name).monitor_function[$idx])
        end
    elseif fp.append
        return :($(fp.name).$name)
    else
        return fp.name
    end
end

function _generate_apply_function(func::Type{CF}, applyto,
                                  args) where {CF <: ContinuationFunction}
    expr = []
    _generate_apply_function!(expr, func, applyto, args)
    return expr
end

function _generate_apply_function!(expr, func::Type{CF}, applyto,
                                   args) where {CF <: ContinuationFunction}
    # ContinuationFunction{Z, M, MNAME, P, PNAME}
    (Z, M, MNAME, P, PNAME) = func.parameters
    for i in Base.OneTo(length(P.parameters))
        name = PNAME.parameters[i]
        new_args = [FuncPar(arg, name, i) for arg in args]
        _generate_apply_function!(expr, P.parameters[i], applyto, new_args)
    end
    if (Z !== Nothing) && ((applyto == :zero) || (applyto == :all))
        apply = Expr(:call)
        for arg in args
            push!(apply.args, get_par(arg, :zero, 0))
        end
        push!(expr, apply)
    end
    if (applyto == :monitor) || (applyto == :all)
        for i in Base.OneTo(length(M.parameters))
            name = MNAME.parameters[i]
            apply = Expr(:call)
            for arg in args
                push!(apply.args, get_par(arg, name, i))
            end
            push!(expr, apply)
        end
    end
    return expr
end

@generated function eval_zero_function!(res, func::ContinuationFunction, u, data, active,
                                        monitor)
    return _generate_eval_zero_function(func)
end

function _generate_eval_zero_function(func::Type{CF}) where {CF <: ContinuationFunction}
    result = quote
        j = 0
    end
    _generate_apply_function!(result.args, func, :zero,
                              [
                                  FuncPar(:func; cont_func = true),
                                  FuncPar(:res; unwrap = true, append = true),
                                  FuncPar(:u; unwrap = true),
                                  FuncPar(:data, unwrap = true),
                              ])
    allmonitor_function = _generate_apply_function(func, :monitor,
                                                   [
                                                       FuncPar(:func; cont_func = true),
                                                       FuncPar(:u; unwrap = true),
                                                       FuncPar(:data, unwrap = true),
                                                       FuncPar(:data, unwrap = true,
                                                               append = true),
                                                   ])
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

@generated function eval_monitor_function!(res, func::ContinuationFunction, u, data)
    return _generate_eval_monitor_function(func)
end

function _generate_eval_monitor_function(func::Type{CF}) where {CF <: ContinuationFunction}
    result = quote end
    allmonitor_function = _generate_apply_function(func, :monitor,
                                                   [
                                                       FuncPar(:func; cont_func = true),
                                                       FuncPar(:u; unwrap = true),
                                                       FuncPar(:data, unwrap = true),
                                                       FuncPar(:data, unwrap = true,
                                                               append = true),
                                                   ])
    for (i, monitor_function) in enumerate(allmonitor_function)
        push!(result.args, :(res[$i] = $monitor_function))
    end
    push!(result.args, :(return res))
    return result
end

"""
    ContinuationData

A wrapper for chart data. Stores current monitor function values and records which
parameters are active.
"""
struct ContinuationData{T, D}
    data::D
    monitor::Vector{T}
    active::Vector{Bool}
end

"""
    ContinuationWrapper

A wrapper for a `ContinuationFunction` that allows it to be used with SciML solvers.
Specifically, it allows `u` and `res` to be regular `Vector`s.

By default `u` and `res` are wrapped in `ComponentVectors`s.
"""
struct ContinuationWrapper{F, T, U, R, M}
    prob::ContinuationProblem
    func::F
    u_wrapper::U
    res_wrapper::R
    monitor_wrapper::M
    monitor_names::Vector{String}
    u0::Vector{T}
    p0::ContinuationData
end

function ContinuationWrapper(prob::ContinuationProblem, pars = nothing)
    return ContinuationWrapper(ContinuationFunction(prob), pars)
end

function ContinuationWrapper(func::ContinuationFunction, pars = nothing)
    prob = func.prob
    (_u0, data) = get_initial(prob)
    u0 = ComponentVector(_u0)
    res_layout = ComponentVector(get_initial_residual_layout(prob))
    monitor = ComponentVector(get_initial_monitor(prob, u0, data))
    active = collect(ComponentVector{Bool}(get_initial_active(prob)))
    monitor_names = monitor_function_name(prob)
    # Activate/deactivate parameters as necessary
    if pars !== nothing
        for par in pars
            if par isa String
                par_name = par
                value = true
            elseif par isa Pair
                par_name = par[1]
                value = par[2]
            else
                throw(ArgumentError("Invalid parameter specification"))
            end
            idx = findfirst(==(par_name), monitor_names)
            if idx === nothing
                throw(ArgumentError("Unknown parameter name"))
            else
                active[idx] = value
            end
        end
    end
    # Extend u0 with monitor function values
    u0_ext = ComponentVector(u0; monitor = monitor[active])
    p0 = ContinuationData(data, collect(monitor), active)
    return ContinuationWrapper(prob, func, Base.Fix2(ComponentVector, getaxes(u0_ext)),
                               Base.Fix2(ComponentVector, getaxes(res_layout)),
                               Base.Fix2(ComponentVector, getaxes(monitor)), monitor_names,
                               collect(u0_ext), p0)
end

function (cw::ContinuationWrapper)(res, u, p::ContinuationData)
    _u = cw.u_wrapper(u)
    _res = cw.res_wrapper(res)
    eval_zero_function!(_res, cw.func, _u, p.data, p.active, p.monitor)
    return res
end

function get_initial(cw::ContinuationWrapper)
    return (cw.u0, cw.p0)
end

get_u_wrapper(cw::ContinuationWrapper) = cw.u_wrapper
get_res_wrapper(cw::ContinuationWrapper) = cw.res_wrapper
get_monitor_wrapper(cw::ContinuationWrapper) = cw.monitor_wrapper
