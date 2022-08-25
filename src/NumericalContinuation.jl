module NumericalContinuation

using DocStringExtensions

const RESERVED_NAMES = Set([:zero])

abstract type AbstractContinuationProblem end

"""
    $TYPEDEF

This is the basic structure required to define a continuation problem. Each
`ContinuationProblem` contains an (optional) zero function, zero or more monitor functions,
and zero or more subproblems. `ContinuationProblem`s are coupled together in a tree
structure to form the overall problem definition.
"""
struct ContinuationProblem <: AbstractContinuationProblem
    zero_function::Any
    monitor_function::Vector{Any}
    monitor_function_names::Vector{Symbol}
    sub_problem::Vector{Any}
    sub_problem_names::Vector{Symbol}
end

function ContinuationProblem(zero_function; monitor = [], sub = [])
    monitor_function = []
    monitor_function_names = Symbol[]
    for m in monitor
        if !(m isa Pair)
            throw(ArgumentError("Monitor functions must be specified in the form `[:name => monitor_function]`"))
        else
            if !(m[1] isa Symbol)
                throw(ArgumentError("Names must be specified as symbols, e.g., `:name`"))
            else
                push!(monitor_function_names, m[1])
                push!(monitor_function, m[2])
            end
        end
    end
    sub_problem = []
    sub_problem_names = Symbol[]
    for p in sub
        if !(p isa Pair)
            throw(ArgumentError("Subproblems must be specified in the form `[:name => sub_problem]`"))
        else
            if !(p[1] isa Symbol)
                throw(ArgumentError("Names must be specified as symbols, e.g., `:name`"))
            else
                push!(sub_problem_names, p[1])
                push!(sub_problem, p[2])
            end
        end
    end
    return ContinuationProblem(zero_function, monitor_function, monitor_function_names,
                               sub_problem, sub_problem_names)
end

"""
    $SIGNATURES

Return the names of the monitor functions in the problem specified and its subproblems.
"""
function monitor_function_names end

function monitor_function_names(prob::ContinuationProblem)
    return _monitor_function_names!(Symbol[], prob, Symbol())
end

function _monitor_function_names!(names::Vector{Symbol}, prob::ContinuationProblem,
                                  prefix)
    for i in eachindex(prob.sub_problem, prob.sub_problem_names)
        name = prob.sub_problem_names[i]
        _monitor_function_names!(names, prob.sub_problem[i], Symbol(prefix, name, :.))
    end
    for i in eachindex(prob.monitor_function, prob.monitor_function_names)
        name = prob.monitor_function_names[i]
        push!(names, Symbol(prefix, name))
    end
    return names
end

"""
    $SIGNATURES

Return the names of the subproblems in the problem specified.
"""
function sub_problem_names end

function sub_problem_names(prob::ContinuationProblem)
    return _sub_problem_names!(Symbol[], prob, Symbol())
end

function _sub_problem_names!(names::Vector{Symbol}, prob::ContinuationProblem,
                             prefix)
    for i in eachindex(prob.sub_problem, prob.sub_problem_names)
        name = prob.sub_problem_names[i]
        _sub_problem_names!(names, prob.sub_problem[i], Symbol(prefix, name, :.))
        push!(names, Symbol(prefix, name))
    end
    return names
end

"""
    $SIGNATURES

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
    $SIGNATURES

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
    $SIGNATURES

Close the `ContinuationProblem` so that structure can no longer be modified (e.g., no new
zero functions / monitor functions can be added). The resulting `ClosedProblem` can be used
efficiently with functions such as `zero_function` and `monitor_function`.
"""
function close_problem end

function close_problem(prob::ContinuationProblem)
    sub_problem = []
    for i in eachindex(prob.sub_problem)
        push!(sub_problem, close_problem(prob.sub_problem[i]))
    end
    return ClosedProblem(prob.zero_function, prob.monitor_function,
                         prob.monitor_function_names, sub_problem, prob.sub_problem_names)
end

"""
    $TYPEDEF

This is a closed (i.e., the structure cannot be modified any further) version of the
`ContinuationProblem` that is used for computational efficiency during continuation.

See [`close_problem`](@ref)
"""
struct ClosedProblem{Z, M, MNAME, P, PNAME} <: AbstractContinuationProblem
    zero_function::Z
    monitor_function::M
    sub_problem::P

    function ClosedProblem(zero_function, monitor_function, monitor_function_names,
                           sub_problem, sub_problem_names)
        if length(monitor_function) != length(monitor_function_names)
            throw(ArgumentError("Each monitor function must have one and only one name"))
        end
        if length(sub_problem) != length(sub_problem_names)
            throw(ArgumentError("Each subproblem must have one and only one name"))
        end
        if !all(isa.(sub_problem, ClosedProblem))
            throw(ArgumentError("All subproblems of a closed problem must also be closed"))
        end
        _check_names([monitor_function_names; sub_problem_names])  # check for duplicates and reserved names
        return new{typeof(zero_function), Tuple{(typeof.(monitor_function))...},
                   Tuple{monitor_function_names...}, Tuple{(typeof.(sub_problem))...},
                   Tuple{sub_problem_names...}}(zero_function, (monitor_function...,),
                                                (sub_problem...,))
    end
end

ClosedProblem(prob::ContinuationProblem) = close_problem(prob)
close_problem(prob::ClosedProblem) = prob

@generated function zero_function(prob::ClosedProblem, u, data)
    result = :([])
    _gen_zero_function!(result.args, prob, :prob, :u, :data)
    return :(reduce(hcat, $result))
end

"""
    $SIGNATURES

INTERNAL ONLY

Generate an expression tree to call `zero_function` on a `ClosedProblem` and its
subproblems, returning a vector. This function should be called from an `@generated`
function. `prob`, `u`, and `data` are the names of the corresponding parameters in the
generated function given as a `Symbol` or `Expr`.
"""
function _gen_zero_function!(result, ::Type{ClosedProblem{Z, M, MNAME, P, PNAME}},
                             prob, u, data) where {Z, M, MNAME, P, PNAME}
    for i in Base.OneTo(length(P.parameters))
        name = PNAME.parameters[i]
        _gen_zero_function!(result, P.parameters[i], :($prob.sub_problem[$i]), :($u.$name),
                            :($data.$name))
    end
    if Z !== Nothing
        push!(result,
              :(zero_function($prob.zero_function, $u.zero, $data.zero;
                              parent = ($u, $data))))
    end
    return result
end

@generated function monitor_function(prob::ClosedProblem, u, data)
    result = :([])
    _gen_monitor_function!(result.args, prob, :prob, :u, :data)
    return result
end

"""
    $SIGNATURES

INTERNAL ONLY

Generate an expression tree to call `monitor_function` on a `ClosedProblem` and its
subproblems, returning a vector. This function should be called from an `@generated`
function. `prob`, `u`, and `data` are the names of the corresponding parameters in the
generated function given as a `Symbol` or `Expr`.
"""
function _gen_monitor_function!(result, ::Type{ClosedProblem{Z, M, MNAME, P, PNAME}},
                                prob, u, data) where {Z, M, MNAME, P, PNAME}
    for i in Base.OneTo(length(P.parameters))
        name = PNAME.parameters[i]
        _gen_monitor_function!(result, P.parameters[i], :($prob.sub_problem[$i]),
                               :($u.$name), :($data.$name))
    end
    for i in Base.OneTo(length(M.parameters))
        name = MNAME.parameters[i]
        push!(result,
              :(monitor_function($prob.monitor_function[$i], $u.zero, $data.$name;
                                 parent = ($u, $data))))
    end
    return result
end

function monitor_function_names(prob::ClosedProblem)
    return _monitor_function_names!(Symbol[], typeof(prob), Symbol())
end

function _monitor_function_names!(names::Vector{Symbol},
                                  ::Type{ClosedProblem{Z, M, MNAME, P, PNAME}},
                                  prefix) where {Z, M, MNAME, P, PNAME}
    # Must match the ordering of _gen_monitor_function!
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

function sub_problem_names(prob::ClosedProblem)
    return _sub_problem_names!(Symbol[], typeof(prob), Symbol())
end

function _sub_problem_names!(names::Vector{Symbol},
                                  ::Type{ClosedProblem{Z, M, MNAME, P, PNAME}},
                                  prefix) where {Z, M, MNAME, P, PNAME}
    # Must match the ordering of _gen_sub_problem!
    for i in Base.OneTo(length(P.parameters))
        name = PNAME.parameters[i]
        _sub_problem_names!(names, P.parameters[i], Symbol(prefix, name, :.))
        push!(names, Symbol(prefix, name))
    end
    return names
end

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
