using QuantumPropagators
using Test
using SafeTestsets

# Note: comment outer @testset to stop after first @safetestset failure
@time @testset verbose = true "QuantumPropagators" begin

    @test v"0.1.0" < QuantumPropagators.VERSION < v"1.0.0"

    println("\n* Doctests:")
    @time @safetestset "Doctests" begin
        using QuantumPropagators, Test, Documenter
        doctest(QuantumPropagators)
    end

    println("\n* Hamiltonian (test_hamiltonian.jl):")
    @time @safetestset "Hamiltonian" begin
        include("test_hamiltonian.jl")
    end

    println("\n* Liouvillian (test_liouvillian.jl):")
    @time @safetestset "Liouvillian" begin
        include("test_liouvillian.jl")
    end

    println("\n* Operator Linalg (test_operator_linalg.jl):")
    @time @safetestset "Operator Linalg" begin
        include("test_operator_linalg.jl")
    end

    println("\n* Shapes (test_shapes.jl):")
    @time @safetestset "Shapes" begin
        include("test_shapes.jl")
    end

    println("\n* Controls (test_controls.jl):")
    @time @safetestset "Controls" begin
        include("test_controls.jl")
    end

    println("\n* Amplitudes (test_amplitudes.jl):")
    @time @safetestset "Amplitudes" begin
        include("test_amplitudes.jl")
    end

    println("\n* Parameterization (test_parameterization.jl):")
    @time @safetestset "Parameterization" begin
        include("test_parameterization.jl")
    end

    println("\n* CRAB Functions (test_crab_functions.jl):")
    @time @safetestset "CRAB Functions" begin
        include("test_crab_functions.jl")
    end

    println("\n* Discretization (test_discretization.jl):")
    @time @safetestset "Discretization" begin
        include("test_discretization.jl")
    end

    println("\n* Spec.Rad. (test_specrad.jl):")
    @time @safetestset "Spec.Rad." begin
        include("test_specrad.jl")
    end

    println("\n* Cheby (test_cheby.jl):")
    @time @safetestset "Cheby" begin
        include("test_cheby.jl")
    end

    println("\n* Newton (test_newton.jl):")
    @time @safetestset "Newton" begin
        include("test_newton.jl")
    end

    println("\n* Exp (test_expprop.jl):")
    @time @safetestset "Exp" begin
        include("test_expprop.jl")
    end

    println("\n* Propagator Interfaces (test_prop_interfaces.jl):")
    @time @safetestset "Propagator Interfaces" begin
        include("test_prop_interfaces.jl")
    end

    println("\n* Propagate (test_propagate.jl):")
    @time @safetestset "Propagate" begin
        include("test_propagate.jl")
    end

    println("\n* Propagate Sequence (test_propagate_sequence.jl):")
    @time @safetestset "Propagate Sequence" begin
        include("test_propagate_sequence.jl")
    end

    println("\n* Time-dependent observables (test_timedependent_observables.jl):")
    @time @safetestset "Time-dependent observables" begin
        include("test_timedependent_observables.jl")
    end

    println("\n* Timings (test_timings.jl):")
    @time @safetestset "Timings" begin
        include("test_timings.jl")
    end

    println("\n* Test non-standard generators (test_nonstandard_generators.jl):")
    @time @safetestset "Non-standard generators" begin
        include("test_nonstandard_generators.jl")
    end

    println("\n* Invalid interfaces (test_invalid_interfaces.jl):")
    @time @safetestset "Invalid interfaces" begin
        include("test_invalid_interfaces.jl")
    end

    print("\n")

end;
nothing
