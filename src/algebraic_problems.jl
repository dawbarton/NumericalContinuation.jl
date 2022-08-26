export algebraic_problem

struct AlgebraicProblem{F}
    f::F
    u0::Any
    p0::Any
end

zero_function(alg::AlgebraicProblem, u, data; parent) = alg.f(u.u, u.p)

function algebraic_problem(f, u0::Union{AbstractVector, NamedTuple},
                           p0::Union{AbstractVector, NamedTuple})
    prob = ContinuationProblem(AlgebraicProblem(f, u0, p0))
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
