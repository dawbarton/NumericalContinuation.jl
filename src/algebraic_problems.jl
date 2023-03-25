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

(alg::AlgebraicProblem)(res, u, data; kwargs...) = alg.f(res, u.zero.u, u.zero.p)

function algebraic_problem(f, u, p = (;), eqns = missing)
    prob = ContinuationProblem(AlgebraicProblem(f, u, p, eqns))
    return add_parameters!(prob, :p, keys(p))
end

@testitem "Algebraic problems" begin
    function alg_problem1(::Type{T}) where {T}
        # Problem 1 (example from COCO book)
        _u0 = T[1.0, 1.0, 0.5]
        prob = algebraic_problem((u, p) -> [u[1]^2 + (u[2] - 1)^2 - 1], _u0)
        add_monitor_function!(prob, :ψ₁, monitor_function((u, p) -> sqrt(u[1]^2 + u[2]^2)))
        add_monitor_function!(prob, :ψ₂, monitor_function((u, p) -> u[1]))
        add_monitor_function!(prob, :ψ₃, monitor_function((u, p) -> u[3] - u[1] + 0.5))
        @test monitor_function_name(prob) == ["ψ₁", "ψ₂", "ψ₃"]
        @test isempty(sub_problem_name(prob))

        # Get initial values
        (u0, data) = NumericalContinuation.get_initial(prob)
        active0 = NumericalContinuation.get_initial_active(prob)
        monitor0 = NumericalContinuation.get_initial_monitor(prob, u0, data)
        res_layout = NumericalContinuation.get_initial_residual_layout(prob)

        # Generate the optimised code
        func = NumericalContinuation.ContinuationFunction(prob)

        # Buffers
        res = ComponentVector{T}(res_layout)
        u = ComponentVector{T}(u0)
        monitor = ComponentVector{T}(monitor0)
        active = ComponentVector{Bool}(active0)
        chart = nothing

        # Test propagated values
        @test u ≈ _u0
        @test monitor ≈ [sqrt(2), 1.0, 0.0]

        # Test with NamedTuple
        monitor1 = zero(monitor)
        NumericalContinuation.eval_monitor_function!(monitor1, func, u0, data)
        @test monitor1 ≈ monitor
        NumericalContinuation.eval_zero_function!(res, func, u0, data, active, monitor)
        @test res ≈ [0, 0, 0, 0]

        # Test with ComponentArray
        NumericalContinuation.eval_monitor_function!(monitor1, func, u, data)
        @test monitor1 ≈ monitor
        NumericalContinuation.eval_zero_function!(res, func, u, data, active, monitor)
        @test res ≈ [0, 0, 0, 0]

        # Test with monitor function active
        u_ext = ComponentVector(u; monitor = [0])
        active[1] = true
        NumericalContinuation.eval_zero_function!(res, func, u_ext, data, active, monitor)
        @test res ≈ [0, sqrt(2), 0, 0]
    end
    alg_problem1(Float64)
    alg_problem1(Float32)
end
