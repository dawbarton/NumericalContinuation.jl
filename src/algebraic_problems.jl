export algebraic_problem

struct AlgebraicProblem{IIP, F}
    f::F
    u0::Any
end

function AlgebraicProblem(f, u0)
    # Try to determine if the problem is in place or not
    for method in methods(f)
        if method.nargs == 4  # (res, u, p) + 1
            return AlgebraicProblem{true, typeof(f)}(f, u0)
        end
    end
    return AlgebraicProblem{false, typeof(f)}(f, u0)
end

(alg::AlgebraicProblem{true})(res, u, data; parent) = alg.f(res, u.u, u.p)
(alg::AlgebraicProblem{false})(res, u, data; parent) = (res .= alg.f(u.u, u.p))

get_initial_data(alg::AlgebraicProblem) = (alg.u0, nothing)

function algebraic_problem(f, u0::Union{AbstractVector, NamedTuple},
                           p0::Union{AbstractVector, NamedTuple})
    prob = ContinuationProblem(AlgebraicProblem(f, ComponentVector((u=u0, p=p0))))
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
