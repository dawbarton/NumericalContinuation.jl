using NumericalContinuation
using Test
using TestItemRunner

@testitem "Aqua.jl tests" begin
    using Aqua: Aqua
    # ComponentArrays has ambiguities...
    Aqua.test_all(NumericalContinuation; ambiguities=false)
end

@run_package_tests
