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

    print("\n* Hamiltonian (test_hamiltonian.jl):")
    @time @safetestset "Hamiltonian" begin
        include("test_hamiltonian.jl")
    end

    print("\n* Liouvillian (test_liouvillian.jl):")
    @time @safetestset "Liouvillian" begin
        include("test_liouvillian.jl")
    end

    print("\n* Operator Linalg (test_operator_linalg.jl):")
    @time @safetestset "Operator Linalg" begin
        include("test_operator_linalg.jl")
    end

    print("\n* Shapes (test_shapes.jl):")
    @time @safetestset "Shapes" begin
        include("test_shapes.jl")
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

    print("\n* Time-dependent observables (test_timedependent_observables.jl):")
    @time @safetestset "Time-dependent observables" begin
        include("test_timedependent_observables.jl")
    end

    print("\n")

end;
