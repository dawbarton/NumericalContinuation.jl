export algebraic_problem

struct AlgebraicProblem
    f::Any
    u0::Any
    vars::Int
    eqns::Int
    function AlgebraicProblem(f, u, p = (;), eqns = missing)
        # If needed, call the function to see how many outputs it returns
        _eqns = ismissing(eqns) ? length(f(u, p)) : eqns
        vars = length(u) + length(p)
        # Empty values for ComponentArrays can only be NamedTuples
        _u = isempty(u) ? (;) : u
        _p = isempty(p) ? (;) : p
        # Ensure that the function has an in-place form
        _f = IIPWrapper(f)
        return new(_f, (u = _u, p = _p), vars, _eqns)
    end
end

(alg::AlgebraicProblem)(res, u, data; kwargs...) = alg.f(res, u.u, u.p)

function algebraic_problem(f, u, p = (;), eqns = missing)
    prob = ContinuationProblem(AlgebraicProblem(f, u, p, eqns))
    return add_parameters!(prob, :p, keys(p))
end

@testitem "Algebraic problems" begin
    function alg_problem1(::Type{T}) where {T}
        # Problem 1 (example from COCO book)
        _u0 = T[1.0, 1.0, 0.5]
        prob = algebraic_problem((u, p) -> [u[1]^2 + (u[2] - 1)^2 - 1], _u0)
        add_monitor_function!(prob, :ψ₁, monitor_function((u, p) -> sqrt(u[1]^2 + u[2]^2); pars = true))
        add_monitor_function!(prob, :ψ₂, monitor_function((u, p) -> u[1]; pars = true))
        add_monitor_function!(prob, :ψ₃, monitor_function((u, p) -> u[3] - u[1] + 0.5; pars = true))
        @test monitor_function_name(prob) == [:ψ₁, :ψ₂, :ψ₃]
        @test isempty(sub_problem_name(prob))
        (u0, data) = NumericalContinuation.get_initial(prob)
        @test (collect(Iterators.flatten(Iterators.flatten(u0))) == _u0)
        res_layout = NumericalContinuation.get_initial_residual_layout(prob)
        @test (length(collect(Iterators.flatten(res_layout))) == 4)

        func = NumericalContinuation.ContinuationFunction(prob)
        monitor = zeros(T, length(monitor_function_name(prob)))
        res = ComponentVector{T}(res_layout)
        # Test with NamedTuple
        NumericalContinuation.eval_monitor_function!(monitor, func, u0, data, nothing)
        @test monitor ≈ [sqrt(2), 1.0, 0.0]
        NumericalContinuation.eval_function!(res, func, u0, data, nothing, Set(Int[]),
                                             monitor)
        @test res ≈ [0, 0, 0, 0]
        # Test with ComponentArray
        u0_cmp = ComponentVector{T}(u0)
        NumericalContinuation.eval_monitor_function!(monitor, func, u0_cmp, data, nothing)
        @test monitor ≈ [sqrt(2), 1.0, 0.0]
        NumericalContinuation.eval_function!(res, func, u0_cmp, data, nothing, Set(Int[]),
                                             monitor)
        @test res ≈ [0, 0, 0, 0]
        # Test with monitor function active
        u0_ext = (; u0..., monitor = [0])
        NumericalContinuation.eval_function!(res, func, u0_ext, data, nothing, Set([1]),
                                             monitor)
        @test res ≈ [0, sqrt(2), 0, 0]
    end
    alg_problem1(Float64)
    alg_problem1(Float32)

    function alg_problem2(::Type{T}) where {T}
        # Problem 1 (example from COCO book) with u[2] := p[1], u[3] := p[2]
        _u0 = T[1.0]
        _p0 = T[1.0, 0.5]
        prob = algebraic_problem((u, p) -> [u[1]^2 + (p[1] - 1)^2 - 1], _u0, _p0)
        add_monitor_function!(prob, :ψ₁, monitor_function((u, p) -> sqrt(u[1]^2 + p[1]^2); pars = true))
        add_monitor_function!(prob, :ψ₂, monitor_function((u, p) -> u[1]; pars = true))
        add_monitor_function!(prob, :ψ₃, monitor_function((u, p) -> p[2] - u[1] + 0.5; pars = true))
        @test monitor_function_name(prob) == [:p1, :p2, :ψ₁, :ψ₂, :ψ₃]
        @test isempty(sub_problem_name(prob))
        (u0, data) = NumericalContinuation.get_initial(prob)
        @test (collect(Iterators.flatten(Iterators.flatten(u0))) == [_u0; _p0])
        res_layout = NumericalContinuation.get_initial_residual_layout(prob)
        @test (length(collect(Iterators.flatten(res_layout))) == 6)

        func = NumericalContinuation.ContinuationFunction(prob)
        monitor = zeros(T, length(monitor_function_name(prob)))
        res = ComponentVector{T}(res_layout)
        # Test with NamedTuple
        NumericalContinuation.eval_monitor_function!(monitor, func, u0, data, nothing)
        @test monitor[3:end] ≈ [sqrt(2), 1.0, 0.0]
        NumericalContinuation.eval_function!(res, func, u0, data, nothing, Set(Int[]),
                                             monitor)
        @test res ≈ [0, 0, 0, 0, 0, 0]
        # Test with ComponentArray
        u0_cmp = ComponentVector{T}(u0)
        NumericalContinuation.eval_monitor_function!(monitor, func, u0_cmp, data, nothing)
        @test monitor[3:end] ≈ [sqrt(2), 1.0, 0.0]
        NumericalContinuation.eval_function!(res, func, u0_cmp, data, nothing, Set(Int[]),
                                             monitor)
        @test res ≈ [0, 0, 0, 0, 0, 0]
        # Test with monitor function active
        u0_ext = (; u0..., monitor = [0])
        NumericalContinuation.eval_function!(res, func, u0_ext, data, nothing, Set([3]),
                                             monitor)
        @test res ≈ [0, 0, 0, sqrt(2), 0, 0]
    end
    alg_problem2(Float64)
    alg_problem2(Float32)
end
