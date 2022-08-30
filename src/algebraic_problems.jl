export algebraic_problem

struct AlgebraicProblem{IIP, F}
    f::F
    u0::Any
    dim::Int
end

function AlgebraicProblem(f, u0; dim::Union{Missing, Integer} = missing)
    # Try to determine if the problem is in place or not
    for method in methods(f)
        if method.nargs == 4  # (res, u, p) + 1
            if dim === missing
                throw(ArgumentError("In place functions must supply the number of dimensions"))
            end
            return AlgebraicProblem{true, typeof(f)}(f, u0, dim)
        end
    end
    if dim === missing
        # Call the function to see how many outputs it returns
        res = f(u0.u, u0.p)
        if res isa AbstractArray
            _f = f
        elseif res isa Number
            _f = (u, p) -> [f(u,p)]
        else
            throw(ArgumentError("Expected the function to return an AbstractArray or a Number"))
        end
        _dim = length(res)
    else
        _dim = dim
    end
    return AlgebraicProblem{false, typeof(_f)}(_f, u0, _dim)
end

(alg::AlgebraicProblem{true})(res, u, data; parent) = alg.f(res, u.u, u.p)
(alg::AlgebraicProblem{false})(res, u, data; parent) = (res .= alg.f(u.u, u.p))

get_initial_data(alg::AlgebraicProblem) = (alg.u0, nothing)
get_residual_vector(alg::AlgebraicProblem, u0, ::Any) = zeros(eltype(u0), alg.dim)

function algebraic_problem(f, u0::Union{AbstractVector, NamedTuple},
                           p0::Union{AbstractVector, NamedTuple} = (;);
                           dim::Union{Missing, Integer} = missing)
    prob = ContinuationProblem(AlgebraicProblem(f, ComponentVector((u = u0, p = p0));
                                                dim = dim))
    # Try to extract names
    if p0 isa NamedTuple
        _p0 = ComponentVector(p0)
    else
        _p0 = p0
    end
    if _p0 isa ComponentVector
        for (i, name) in enumerate(labels(_p0))
            add_parameter!(prob, Symbol(name), @optic _.p[i])
        end
    else
        add_parameter!(prob, 1:length(p0))
    end
    return prob
end

# prob = algebraic_problem((u, p) -> [u.x^2 + u.y^2 - p.r], (x=1.0, y=0.0), (r=1.0,))

function test1()
    prob = algebraic_problem((u, p) -> u[1]^2 + (u[2] - 1)^2 - 1, [1.0, 1.0, 0.5])
    add_monitor_function!(prob, :ψ₁, monitor_function(u -> sqrt(u[1]^2 + u[2]^2)))
    add_monitor_function!(prob, :ψ₂, monitor_function(u -> u[1]))
    add_monitor_function!(prob, :ψ₃, monitor_function(u -> u[3] - u[1] + 0.5))
    return prob
end
