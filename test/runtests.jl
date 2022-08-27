using QuantumPropagators
using Test
using SafeTestsets

# Note: comment outer @testset to stop after first @safetestset failure
@time @testset verbose = true "QuantumPropagators" begin

    print("\n* Doctests:")
    @time @safetestset "Doctests" begin
        using QuantumPropagators, Test, Documenter
        doctest(QuantumPropagators)
    end

    print("\n* Controls (test_controls.jl):")
    @time @safetestset "Controls" begin
        include("test_controls.jl")
    end

    print("\n* Discretization (test_discretization.jl):")
    @time @safetestset "Discretization" begin
        include("test_discretization.jl")
    end

    print("\n* Spec.Rad. (test_specrad.jl):")
    @time @safetestset "Spec.Rad." begin
        include("test_specrad.jl")
    end

    print("\n* Cheby (test_cheby.jl):")
    @time @safetestset "Cheby" begin
        include("test_cheby.jl")
    end

    print("\n* Newton (test_newton.jl):")
    @time @safetestset "Newton" begin
        include("test_newton.jl")
    end

    print("\n* Exp (test_expprop.jl):")
    @time @safetestset "Exp" begin
        include("test_expprop.jl")
    end

    print("\n* Propagate (test_propagate.jl):")
    @time @safetestset "Propagate" begin
        include("test_propagate.jl")
    end

    print("\n")

end;
