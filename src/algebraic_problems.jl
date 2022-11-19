export algebraic_problem

struct AlgebraicProblem{F}
    f::F
    u0::Any
    dim::Int
    function AlgebraicProblem(f, u0, p0=(;), dim=missing)
        if ismissing(dim)
            # Call the function to see how many outputs it returns
            res = f(u0, p0)
            if !(res isa AbstractArray)
                throw(ArgumentError("Expected the function to return an AbstractArray"))
            end
            dim = length(res)
        end
        # Empty values for ComponentArrays can only be NamedTuples
        _u0 = isempty(u0) ? (;) : u0
        _p0 = isempty(p0) ? (;) : p0
        # Ensure that the function has an in-place form
        _f = IIPWrapper(f)
        return new{typeof(_f)}(_f, ComponentVector((u = _u0, p = _p0)), dim)
    end
end

(alg::AlgebraicProblem)(res, u, data; parent) = alg.f(res, u.u, u.p)

get_initial_data(alg::AlgebraicProblem) = (alg.u0, nothing)
get_residual_vector(alg::AlgebraicProblem, u0, ::Any) = zeros(eltype(u0), alg.dim)

function algebraic_problem(f, u0, p0=(;), dim=missing)
    prob = ContinuationProblem(AlgebraicProblem(f, u0, p0, dim))
    return add_parameter_p0!(prob, p0)
end

@testitem "Algebraic problems" begin
    function alg_problem1(::Type{T}) where T
        # Problem 1 (example from COCO book)
        _u0 = T[1.0, 1.0, 0.5]
        prob = algebraic_problem((u, p) -> [u[1]^2 + (u[2] - 1)^2 - 1], _u0)
        add_monitor_function!(prob, :ψ₁, monitor_function(u -> sqrt(u[1]^2 + u[2]^2)))
        add_monitor_function!(prob, :ψ₂, monitor_function(u -> u[1]))
        add_monitor_function!(prob, :ψ₃, monitor_function(u -> u[3] - u[1] + 0.5))
        @test monitor_function_name(prob) == [:ψ₁, :ψ₂, :ψ₃]
        @test isempty(sub_problem_name(prob))
        (u0, data) = NumericalContinuation.get_initial_data(prob)
        @test (collect(u0) == _u0) && (eltype(u0) == eltype(_u0))
        res = NumericalContinuation.get_residual_vector(prob, u0, data)
        @test (length(res) == 4) && (eltype(res) == eltype(_u0))
        func = NumericalContinuation.ContinuationFunction(prob)
        monitor = zeros(length(monitor_function_name(prob)))
        NumericalContinuation.eval_monitor_function!(monitor, func, u0, data)
        @test monitor ≈ [sqrt(2), 1.0, 0.0]
        res .= one(eltype(res))
        NumericalContinuation.eval_function!(res, func, u0, data, Set(Int[]), monitor)
        @test res ≈ [0, 0, 0, 0]
        u0_ext = [u0; ComponentVector(monitor=[0])]
        NumericalContinuation.eval_function!(res, func, u0_ext, data, Set([1]), monitor)
        @test res ≈ [0, sqrt(2), 0, 0]
    end
    alg_problem1(Float64)
    alg_problem1(Float32)
end
