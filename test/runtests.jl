using QuantumPropagators
using Test
using SafeTestsets

# Note: comment outer @testset to stop after first @safetestset failure
@time @testset verbose=true "QuantumPropagators" begin

    print("\n**** Chebychev propagator")
    @time @safetestset "Chebychev propagator" begin include("test_cheby.jl") end

    print("\n**** Newton propagator")
    @time @safetestset "Newton propagator" begin include("test_newton.jl") end

    print("\n**** Exp propagator")
    @time @safetestset "Exp propagator" begin include("test_expprop.jl") end

    print("\n**** Total\n")

end
