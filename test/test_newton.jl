using Test
using LinearAlgebra
using QuantumPropagators

include("utils.jl")

@testset "random Hermitian" begin

    SERIALIZATION = 0
    SERIALIZE_NUMPY = true
    # A normal test run is with SERIALIZATION = 0. This uses a random
    # Hamiltonian and initial state
    # Set NORMALIZATION = 1 to write the random Hamiltonian and initial state
    # to disk, and then NORMALIZATION = 2 to load it from disk (ensuring a
    # reproducible test run.
    # In combination with NORMALIZATION = 1, if SERIALIZE_NUMPY is true, the
    # Hamiltonian and initial state will also be dumped into a numpy (.npy)
    # file, which then can be read into the test_newton.py script.

    precision = 1e-10

    # Input
    N = 1000
    ρ = 10  # approximate spectral radius

    dt = 0.5

    if SERIALIZATION < 2
        H = random_hermitian_matrix(N, ρ)
        Ψ₀ = random_state_vector(N)
    end

    if SERIALIZATION == 1
        serialize("test/H.dump", H)
        serialize("test/psi0.dump", Ψ₀)
        if SERIALIZE_NUMPY
            using PyCall
            np = pyimport("numpy")
            # write numpy dump files psi0.npy, H.npy;
            # can be loaded with `np.load`
            np.save("test/psi0", Ψ₀)
            np.save("test/H", H)
        end
    end
    if SERIALIZATION == 2
        H = deserialize("test/H.dump")
        Ψ₀ = deserialize("test/psi0.dump")
    end

    @test norm(Ψ₀) ≈ 1

    # Expected
    U = exp(-im * H * dt)
    @test norm(U * U'  - one(U)) < precision  # U is unitary
    Ψ_out_expected = U * Ψ₀
    @test norm(Ψ_out_expected) ≈ 1

    Ψ = copy(Ψ₀)
    wrk = NewtonWrk(Ψ₀; m_max=5)
    newton!(Ψ, H, dt, wrk; max_restarts=200)
    Ψ_out = copy(Ψ)
    @test norm(Ψ_out) ≈ 1

    # Comparison
    @test norm(Ψ_out - Ψ_out_expected) < precision

end


@testset "random non-Hermitian" begin

    SERIALIZATION = 0
    SERIALIZE_NUMPY = true
    # A normal test run is with SERIALIZATION = 0. This uses a random
    # Hamiltonian and initial state
    # Set NORMALIZATION = 1 to write the random Hamiltonian and initial state
    # to disk, and then NORMALIZATION = 2 to load it from disk (ensuring a
    # reproducible test run.
    # In combination with NORMALIZATION = 1, if SERIALIZE_NUMPY is true, the
    # Hamiltonian and initial state will also be dumped into a numpy (.npy)
    # file, which then can be read into the test_newton.py script.

    precision = 1e-10

    # Input
    N = 1000
    ρ = 10  # approximate spectral radius

    dt = 0.5

    if SERIALIZATION < 2
        H = random_complex_matrix(N, ρ)
        Ψ₀ = random_state_vector(N)
    end

    if SERIALIZATION == 1
        serialize("test/H.dump", H)
        serialize("test/psi0.dump", Ψ₀)
        if SERIALIZE_NUMPY
            using PyCall
            np = pyimport("numpy")
            # write numpy dump files psi0.npy, H.npy;
            # can be loaded with `np.load`
            np.save("test/psi0", Ψ₀)
            np.save("test/H", H)
        end
    end
    if SERIALIZATION == 2
        H = deserialize("test/H.dump")
        Ψ₀ = deserialize("test/psi0.dump")
    end

    @test norm(Ψ₀) ≈ 1

    # Expected
    U = exp(-im * H * dt)
    Ψ_out_expected = U * Ψ₀

    Ψ = copy(Ψ₀)
    wrk = NewtonWrk(Ψ₀; m_max=50)
    newton!(Ψ, H, dt, wrk; max_restarts=200)
    Ψ_out = copy(Ψ)

    # Comparison
    @test norm(Ψ_out - Ψ_out_expected) < precision

end


@testset "random sparse Liouvillian" begin

    SERIALIZATION = 0
    # A normal test run is with SERIALIZATION = 0. This uses a random
    # Hamiltonian and initial state
    # Set NORMALIZATION = 1 to write the random Hamiltonian and initial state
    # to disk, and then NORMALIZATION = 2 to load it from disk (ensuring a
    # reproducible test run.

    precision = 1e-10

    # Input
    N = 32
    ρ = 10  # approximate spectral radius
    sparsity = 0.5

    dt = 0.5

    if SERIALIZATION < 2
        L = random_complex_sparse_matrix(N^2, ρ, sparsity)
        Ψ₀ = random_state_vector(N)
        ρ₀ = reshape(Ψ₀ * Ψ₀', :)
    end

    if SERIALIZATION == 1
        serialize("test/H.dump", L)
        serialize("test/rho0.dump", ρ₀)
    end
    if SERIALIZATION == 2
        L = deserialize("test/H.dump")
        ρ₀ = deserialize("test/psi0.dump")
    end

    @test tr(reshape(ρ₀, (N, N))) ≈ 1

    # Expected
    U = exp(Array(-im * L * dt))
    ρ_out_expected = U * ρ₀

    ρ = copy(ρ₀)
    wrk = NewtonWrk(ρ₀; m_max=50)
    newton!(ρ, L, dt, wrk; max_restarts=20)
    ρ_out = copy(ρ)

    # Comparison
    @test norm(ρ_out - ρ_out_expected) < precision

end
