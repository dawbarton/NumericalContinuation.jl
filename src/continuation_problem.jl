# Types
export ContinuationProblem, ContinuationFunction
# Functions
export add_monitor_function!, add_sub_problem!, monitor_function_name, sub_problem_name,
       add_parameter!, monitor_function
# Re-exports
export @optic, ComponentVector

"""
    ContinuationProblem

This is the basic structure required to define a continuation problem. Each
`ContinuationProblem` contains an (optional) zero function, zero or more monitor functions,
and zero or more subproblems. `ContinuationProblem`s are coupled together in a tree
structure to form the overall problem definition.
"""
struct ContinuationProblem
    zero_function!::Any
    monitor_function::Vector{Any}
    monitor_function_name::Vector{Symbol}
    sub_problem::Vector{Any}
    sub_problem_name::Vector{Symbol}
end

"""
    ContinuationProblem(zero_function!; [monitor], [subprob])

Create a `ContinuationProblem` from a zero function, optionally adding in monitor functions
and subproblems.
"""
function ContinuationProblem(zero_function!; monitor = [], subprob = [])
    prob = ContinuationProblem(zero_function!, [], Symbol[], [], Symbol[])
    isempty(subprob) || add_sub_problem!(prob, subprob)
    isempty(monitor) || add_monitor_function!(prob, monitor)
    return prob
end

"""
    add_monitor_function!(prob, name, mfunc)
    add_monitor_function!(prob, named_mfunc)
    add_monitor_function!(prob, named_mfuncs)

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
    add_sub_problem!(prob, name, sub_prob)
    add_sub_problem!(prob, named_sub_prob)
    add_sub_problem!(prob, named_sub_probs)

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
    monitor_function_name(prob)

Return the names of the monitor functions in the problem specified and its subproblems as a
`Vector{String}`.
"""
function monitor_function_name end

function monitor_function_name(prob::ContinuationProblem)
    return _monitor_function_names!(String[], prob, "")
end

function _monitor_function_names!(names::Vector{String}, prob::ContinuationProblem,
                                  prefix)
    # Must match the ordering of _generate_apply_function
    for i in eachindex(prob.sub_problem_name)
        name = prob.sub_problem_name[i]
        _monitor_function_names!(names, prob.sub_problem[i], "$(prefix)$(name).")
    end
    for i in eachindex(prob.monitor_function_name)
        name = prob.monitor_function_name[i]
        push!(names, "$(prefix)$(name)")
    end
    return names
end

"""
    sub_problem_name(prob)

Return the names of the subproblems in the problem specified as a `Vector{String}`.
"""
function sub_problem_name end

function sub_problem_name(prob::ContinuationProblem)
    return _sub_problem_names!(String[], prob, "")
end

function _sub_problem_names!(names::Vector{String}, prob::ContinuationProblem,
                             prefix)
    for i in eachindex(prob.sub_problem_name)
        name = prob.sub_problem_name[i]
        _sub_problem_names!(names, prob.sub_problem[i], "$(prefix)$(name).")
        push!(names, "$(prefix)$(name)")
    end
    return names
end

"""
    get_vars(prob)

Return the number of variables defined in a problem. Does not include monitor variables.
"""
function get_vars(prob::ContinuationProblem)
    return reduce(sum ∘ get_vars, prob.sub_problem; init = 0) +
           get_vars(prob.zero_function!)
end

get_vars(zero) = zero.vars  # default fallback for zero functions

"""
    get_eqns(prob)

Return the number of equations defined in a problem, including monitor functions.
"""
function get_eqns(prob::ContinuationProblem)
    return reduce(sum ∘ get_eqns, prob.sub_problem; init = 0) +
           get_eqns(prob.zero_function!) +
           length(prob.monitor_function)
end

get_eqns(zero) = zero.eqns  # default fallback for zero functions

"""
    get_eltype(prob)

Return the eltype of the initial state of a problem.
"""
function get_eltype(prob::ContinuationProblem)
    # Could use `reduce(promote_type, map(get_eltype, prob.sub_problem); init=Union{})` but it allocates
    T = Union{}
    # Sub-problems
    for sub_prob in prob.sub_problem
        T = promote_type(get_eltype(sub_prob), T)
    end
    # Zero function
    return promote_type(get_eltype(prob, prob.zero_function!), T)
end

get_eltype(::ContinuationProblem, zero) = get_eltype(zero.u0)  # default fallback for zero function

get_eltype(arr) = eltype(arr)

function get_eltype(tpl::Union{Tuple, NamedTuple})
    reduce(promote_type, map(get_eltype, tpl); init = Union{})
end

"""
    get_initial(prob)

Return a tuple of the initial values of the state vector, and chart data, as (potentially
nested) `NamedTuple`s.
"""
function get_initial(prob::ContinuationProblem)
    u = Pair{Symbol, Any}[]
    data = Pair{Symbol, Any}[]
    # Sub-problems
    for (name, sub_prob) in zip(prob.sub_problem_name, prob.sub_problem)
        (_u, _data) = get_initial(sub_prob)
        push!(u, name => _u)
        push!(data, name => _data)
    end
    # Zero function
    (_u, _data) = get_initial(prob, prob.zero_function!)
    push!(u, :zero => _u)
    push!(data, :zero => _data)
    # Monitor functions
    for (name, monitor_function) in zip(prob.monitor_function_name, prob.monitor_function)
        push!(data, name => get_initial(prob, monitor_function, name))
    end
    return (NamedTuple(u), NamedTuple(data))
end

get_initial(::ContinuationProblem, zero) = (zero.u0, ())  # default fallback for zero function
get_initial(::ContinuationProblem, monitor, name) = ()  # default fallback for monitor functions

"""
    get_initial_residual_layout(prob)

Return a (potentially nested) `NamedTuple` of zeros in the required shape for the residual.
"""
function get_initial_residual_layout(prob::ContinuationProblem)
    (res, monitor) = _get_initial_residual_layout(prob)
    return merge(res, (monitor = falses(monitor),))
end

function _get_initial_residual_layout(prob::ContinuationProblem)
    res = []
    monitor = 0
    for (name, sub_prob) in zip(prob.sub_problem_name, prob.sub_problem)
        (_res, _monitor) = _get_initial_residual_layout(sub_prob)
        push!(res, name => _res)
        monitor += _monitor
    end
    push!(res, :zero => falses(get_eqns(prob.zero_function!)))
    monitor += length(prob.monitor_function)
    return (NamedTuple(res), monitor)
end

"""
    get_initial_monitor(prob, u, data)

Return a tuple of the initial values of the monitor functions as (potentially nested)
`NamedTuple`s.
"""
function get_initial_monitor(prob::ContinuationProblem, u, data)
    monitor = Pair{Symbol, Any}[]
    # Sub-problems
    for (name, sub_prob) in zip(prob.sub_problem_name, prob.sub_problem)
        _monitor = get_initial_monitor(sub_prob, getproperty(u, name),
                                       getproperty(data, name))
        push!(monitor, name => _monitor)
    end
    # Monitor functions
    for (name, monitor_function) in zip(prob.monitor_function_name, prob.monitor_function)
        _monitor = get_initial_monitor(prob, monitor_function, name, u, data)
        push!(monitor, name => _monitor)
    end
    return NamedTuple(monitor)
end

function get_initial_monitor(::ContinuationProblem, mfunc, name, u, data)
    # Default fallback for monitor functions that follow the MonitorFunction interface
    if mfunc.initial_value !== nothing
        value = mfunc.initial_value
    else
        # Call the monitor function to get its initial value
        value = mfunc(u, data, getproperty(data, name))
    end
    return value
end

"""
    get_initial_active(prob, u, data)

Return a tuple of the initial values of the monitor functions, and whether they are active
or not, as (potentially nested) `NamedTuple`s.
"""
function get_initial_active(prob::ContinuationProblem)
    active = Pair{Symbol, Any}[]
    # Sub-problems
    for (name, sub_prob) in zip(prob.sub_problem_name, prob.sub_problem)
        _active = get_initial_active(sub_prob)
        push!(active, name => _active)
    end
    # Monitor functions
    for (name, monitor_function) in zip(prob.monitor_function_name, prob.monitor_function)
        _active = get_initial_active(prob, monitor_function)
        push!(active, name => _active)
    end
    return NamedTuple(active)
end

# Default fallback for monitor functions that follow the MonitorFunction interface
get_initial_active(::ContinuationProblem, mfunc) = mfunc.initial_active

struct MonitorFunction{D, F}
    f::F
    initial_value::Any
    initial_active::Bool
    function MonitorFunction(f, value, active, data = false)
        return new{data, typeof(f)}(f, value, active)
    end
end

(monitor::MonitorFunction{false})(u, data, mdata; kwargs...) = monitor.f(u)

function (monitor::MonitorFunction{true})(u, data, mdata; kwargs...)
    return monitor.f(u, data, mdata)
end

"""
    monitor_function(f; value = nothing, active = false, data = false)

Wrap a function of the form `f(u)` in a monitor function-compatible form. The initial value
of the monitor function may be provided in `value` and whether the monitor function is
active is determined by `active`.

If `value` is not provided, the monitor function will be called with the initial state to
determine the initial value.

If chart data is required, `data` should be set to `true`. In this case, the monitor
function will be called with the form `f(u, data, mdata)`, where `data` is all the problem
data and `mdata` is the monitor function data.
"""
function monitor_function(f; value = nothing, active = false, data = false)
    return MonitorFunction(f, value, active, data)
end

"""
    add_parameter!(prob, name, func; active = false, value = nothing)

Add a continuation parameter `name` to a continuation problem. Parameters should be
specified as function that takes the entire state vector and returns the parameter.

Parameters can be marked as active or inactive. If a parameter is inactive, it will not be
allowed to change from the initial value. If a parameter is active, it will be allowed to
change from the initial value. The initial value of the parameter can be specified in
`value`.

It is assumed that the parameters are part of the state vector from the zero function and so
it is not possible to access variables from the sub-problems. If this is required, use
`monitor_function` directly.

# Examples

```julia
add_parameter!(prob, :α, @optic(_[7]))  # parameter α corresponds to u[7]
add_parameter!(prob, :p, @optic(_.p[1]))  # parameter p corresponds to u.p[1]
```
"""
function add_parameter! end

function add_parameter!(prob::ContinuationProblem, name::Symbol, lens; active = false,
                        value = nothing)
    add_monitor_function!(prob, name,
                          monitor_function(lens ∘ (@optic _.zero); active = active,
                                           value = value))
end

function add_parameter!(prob::ContinuationProblem, name::Symbol, idx::Symbol; kwargs...)
    add_parameter!(prob, Symbol(idx), opcompose(PropertyLens(name), PropertyLens(idx));
                   kwargs...)
end

function add_parameter!(prob::ContinuationProblem, name::Symbol, idx::Integer; kwargs...)
    add_parameter!(prob, Symbol(idx), opcompose(PropertyLens(name), IndexLens(idx));
                   kwargs...)
end

"""
    add_parameters!(prob, name, indices; active = false)

A convenience wrapper for `add_parameter!` that adds multiple continuation parameters at
once. It automatically constructs the parameter mapping based on the name and indices
provided.

It is assumed that the parameters are part of the state vector from the zero function and so
it is not possible to access variables from the sub-problems. If this is required, use
`monitor_function` directly.

# Examples

```julia
add_parameters!(prob, :p, 1:2)  # adds p1 = u.p[1] and p2 = u.p[2]
add_parameters!(prob, :p, [:a, :b])  # adds a = u.p.a and b = u.p.b
```
"""
function add_parameters!(prob::ContinuationProblem, name::Symbol, indices; active = false)
    for idx in indices
        add_parameter!(prob, name, idx; active = active)
    end
    return prob
end

"""
    IIPWrapper

A function wrapper to ensure that out-of-place functions are wrapped as is-in-place
functions.
"""
struct IIPWrapper{IIP, F}
    f::F
end

IIPWrapper(f::IIPWrapper) = f

"""
    IIPWrapper(f, [n=2])

Provide a wrapper to generate in-place versions of a function. `n` is the number of expected
arguments for a non-inplace version. (E.g., `f(u, p)` gives `n=2`.)
"""
function IIPWrapper(f, n = 2)
    # Try to determine if the problem is in place or not
    for method in methods(f)
        if method.nargs == n + 2
            return IIPWrapper{true, typeof(f)}(f)
        end
    end
    return IIPWrapper{false, typeof(f)}(f)
end

"""
    is_iip(func::IIPWrapper)

Returns true if the wrapper function is-in-place or false if it is out-of-place.
"""
is_iip(::IIPWrapper{IIP}) where {IIP} = IIP

(func::IIPWrapper{false})(res, args...; kwargs...) = (res .= func.f(args...; kwargs...))
(func::IIPWrapper{true})(res, args...; kwargs...) = func.f(res, args...; kwargs...)
