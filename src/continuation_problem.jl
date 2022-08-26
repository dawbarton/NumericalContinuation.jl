# Types
export ContinuationProblem, ClosedProblem
# Functions
export add_monitor_function!, add_sub_problem!, monitor_function_name, sub_problem_name,
       add_parameter!
# Re-exports
export @optic, ComponentVector

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
    monitor_function_name::Vector{Symbol}
    sub_problem::Vector{Any}
    sub_problem_name::Vector{Symbol}
end

function ContinuationProblem(zero_function; monitor = [], sub = [])
    prob = ContinuationProblem(zero_function, [], Symbol[], [], Symbol[])
    isempty(sub) || add_sub_problem!(prob, sub)
    isempty(monitor) || add_monitor_function!(prob, monitor)
    return prob
end

"""
    $SIGNATURES

Add named monitor function(s) to a continuation problem.

# Examples

```julia
add_monitor_function!(prob, :some_name, my_monitor_func)
add_monitor_function!(prob, :another_name=>my_monitor_func)
add_monitor_function!(prob, [:mfunc1=>mfunc1, :mfunc2=>mfunc2])
```
"""
function add_monitor_function! end

function add_monitor_function!(prob::ContinuationProblem, name::Symbol, mfunc)
    push!(prob.monitor_function, mfunc)
    push!(prob.monitor_function_name, name)
    return prob
end

function add_monitor_function!(prob::ContinuationProblem, named_mfunc::Pair{Symbol})
    return add_monitor_function!(prob, named_mfunc[1], named_mfunc[2])
end

function add_monitor_function!(prob::ContinuationProblem, named_mfuncs)
    for named_mfunc in named_mfuncs
        if !(named_mfunc isa Pair{Symbol})
            throw(ArgumentError("Monitor functions must be specified in the form `[:name => monitor_function]`"))
        else
            add_monitor_function!(prob, named_mfunc)
        end
    end
    return prob
end

"""
    $SIGNATURES

Add named subproblems function(s) to a continuation problem.

# Examples

```julia
add_sub_problem!(prob, :some_name, my_sub_problem)
add_sub_problem!(prob, :another_name=>my_sub_problem)
add_sub_problem!(prob, [:subprob1=>subprob1, :subprob2=>subprob2])
```
"""
function add_sub_problem! end

function add_sub_problem!(prob::ContinuationProblem, name::Symbol, sub_prob)
    push!(prob.sub_problem, sub_prob)
    push!(prob.sub_problem_name, name)
    return prob
end

function add_sub_problem!(prob::ContinuationProblem, named_sub_prob::Pair{Symbol})
    return add_sub_problem!(prob, named_sub_prob[1], named_sub_prob[2])
end

function add_sub_problem!(prob::ContinuationProblem, named_sub_probs)
    for named_sub_prob in named_sub_probs
        if !(named_sub_prob isa Pair{Symbol})
            throw(ArgumentError("Subproblems must be specified in the form `[:name => sub_problem]`"))
        else
            add_sub_problem!(prob, named_sub_prob)
        end
    end
    return prob
end

"""
    $SIGNATURES

Return the names of the monitor functions in the problem specified and its subproblems.
"""
function monitor_function_name end

function monitor_function_name(prob::AbstractContinuationProblem)
    return _monitor_function_names!(Symbol[], prob, Symbol())
end

function _monitor_function_names!(names::Vector{Symbol}, prob::AbstractContinuationProblem,
                                  prefix)
    # Must match the ordering of _gen_monitor_function!
    for i in eachindex(prob.sub_problem, prob.sub_problem_name)
        name = prob.sub_problem_name[i]
        _monitor_function_names!(names, prob.sub_problem[i], Symbol(prefix, name, :.))
    end
    for i in eachindex(prob.monitor_function, prob.monitor_function_name)
        name = prob.monitor_function_name[i]
        push!(names, Symbol(prefix, name))
    end
    return names
end

"""
    $SIGNATURES

Return the names of the subproblems in the problem specified.
"""
function sub_problem_name end

function sub_problem_name(prob::AbstractContinuationProblem)
    return _sub_problem_names!(Symbol[], prob, Symbol())
end

function _sub_problem_names!(names::Vector{Symbol}, prob::AbstractContinuationProblem,
                             prefix)
    for i in eachindex(prob.sub_problem, prob.sub_problem_name)
        name = prob.sub_problem_name[i]
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
                         prob.monitor_function_name, sub_problem, prob.sub_problem_name)
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
    monitor_function_name::Vector{Symbol}
    sub_problem::P
    sub_problem_name::Vector{Symbol}

    function ClosedProblem(zero_function, monitor_function, monitor_function_name,
                           sub_problem, sub_problem_name)
        if length(monitor_function) != length(monitor_function_name)
            throw(ArgumentError("Each monitor function must have one and only one name"))
        end
        if length(sub_problem) != length(sub_problem_name)
            throw(ArgumentError("Each subproblem must have one and only one name"))
        end
        if !all(isa.(sub_problem, ClosedProblem))
            throw(ArgumentError("All subproblems of a closed problem must also be closed"))
        end
        _check_names([monitor_function_name; sub_problem_name])  # check for duplicates and reserved names
        return new{typeof(zero_function), Tuple{(typeof.(monitor_function))...},
                   Tuple{monitor_function_name...}, Tuple{(typeof.(sub_problem))...},
                   Tuple{sub_problem_name...}}(zero_function, (monitor_function...,),
                                                copy(monitor_function_name),
                                                (sub_problem...,), copy(sub_problem_name))
    end
end

ClosedProblem(prob::ContinuationProblem) = close_problem(prob)
close_problem(prob::ClosedProblem) = prob

@generated function eval_zero_function(prob::ClosedProblem, u, data)
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
        push!(result, :($prob.zero_function($u.zero, $data.zero; parent = ($u, $data))))
    end
    return result
end

@generated function eval_monitor_function(prob::ClosedProblem, u, data)
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
              :($prob.monitor_function[$i]($u.zero, $data.$name; parent = ($u, $data))))
    end
    return result
end

function get_initial_data(prob::AbstractContinuationProblem)
    (u0, data) = _get_initial_data(prob)
    return (ComponentVector(u0), data)
end

function _get_initial_data(prob::AbstractContinuationProblem)
    u0 = []
    data = []
    for i in eachindex(prob.sub_problem, prob.sub_problem_name)
        name = prob.sub_problem_name[i]
        (_u0, _data) = _get_initial_data(prob.sub_problem[i])
        push!(u0, name => _u0)
        push!(data, name => _data)
    end
    (_u0, _data) = get_initial_data(prob.zero_function)
    push!(u0, :zero => _u0)
    push!(data, :zero => _data)
    for i in eachindex(prob.monitor_function, prob.monitor_function_name)
        name = prob.monitor_function_name[i]
        # TODO: decide if monitor functions can introduce new state variables
        _data = get_initial_data(prob.monitor_function[i])
        push!(data, name => _data)
    end
    return (NamedTuple(u0), NamedTuple(data))
end

get_initial_data(monitor_function) = nothing  # fall back for monitor functions

"""
    ContinuationParameter(name, lens)

Create a continuation parameter for use as a monitor function. The reference to a parameter
is stored as a lens (see Accessors.jl) meaning that arbitrary references are possible.

# Example

```julia
ContinuationParameter(:mypar, @optic)
```

Also see: [`add_parameter!`](@ref)
"""
struct ContinuationParameter{L}
    name::Symbol
    lens::L
end

function (cpar::ContinuationParameter)(u, data; parent)
    return cpar.lens(u)
end

"""
    $(SIGNATURES)

Add a continuation parameter to a continuation problem. Parameters should be specified as
`:name=>index` pairs; if only names are provided, the index is taken from the position in
the vector.
"""
function add_parameter! end

function add_parameter!(prob::ContinuationProblem, name::Symbol, lens)
    add_monitor_function!(prob, name, ContinuationParameter(name, lens))
end

function add_parameter!(prob::ContinuationProblem, name::Symbol, idx::Integer)
    add_parameter!(prob, name, @optic _[idx])
end

function add_parameter!(prob::ContinuationProblem, named_pair::Pair{Symbol})
    add_parameter!(prob, named_pair[1], named_pair[2])
end

function add_parameter!(prob::ContinuationProblem, names)
    for (i, name) in enumerate(names)
        if name isa Symbol
            add_parameter!(prob, name, i)
        elseif name isa Integer
            add_parameter!(prob, Symbol(:p, name), name)
        elseif name isa Pair{Symbol}
            add_parameter!(prob, name)
        else
            throw(ArgumentError("Parameters should be specified in the form `[:name=>index]`"))
        end
    end
    return prob
end
