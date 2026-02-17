using Test
using LinearAlgebra
using SparseArrays
using QuantumPropagators
using QuantumPropagators.Generators: liouvillian
using ExponentialUtilities
using StableRNGs: StableRNG


@testset "Forward propagation (Hermitian)" begin
    # Simple TLS Rabi cycling: compare with exact matrix exponential
    Ψ₀ = ComplexF64[1, 0]
    Ĥ = ComplexF64[0 0.5; 0.5 0]
    tlist = collect(range(0, 1.5π, length = 101))
    generator = (Ĥ,)

    Ψ_out = propagate(Ψ₀, generator, tlist; method = ExponentialUtilities)
    Ψ_expected = ComplexF64[-1/√2, -1im/√2]

    @test norm(Ψ_out - Ψ_expected) < 1e-9
end


@testset "Forward propagation (not in-place)" begin
    Ψ₀ = ComplexF64[1, 0]
    Ĥ = ComplexF64[0 0.5; 0.5 0]
    tlist = collect(range(0, 1.5π, length = 101))
    generator = (Ĥ,)

    Ψ_out = propagate(Ψ₀, generator, tlist; method = ExponentialUtilities, inplace = false)
    Ψ_expected = ComplexF64[-1/√2, -1im/√2]

    @test norm(Ψ_out - Ψ_expected) < 1e-9
end


@testset "Backward propagation" begin
    Ψ₀ = ComplexF64[1, 0]
    Ĥ = ComplexF64[0 0.5; 0.5 0]
    tlist = collect(range(0, 1.5π, length = 101))
    generator = (Ĥ,)

    Ψ_fw = propagate(Ψ₀, generator, tlist; method = ExponentialUtilities)
    Ψ_bw = propagate(Ψ_fw, generator, tlist; method = ExponentialUtilities, backward = true)

    @test norm(Ψ_bw - Ψ₀) < 1e-9
end


@testset "Liouvillian propagation" begin
    rng = StableRNG(20260216)
    # 2-level system with decay (Liouvillian)
    σ_x = ComplexF64[0 1; 1 0]
    σ_z = ComplexF64[1 0; 0 -1]
    γ = 0.1  # decay rate
    c_op = √γ * ComplexF64[0 1; 0 0]  # lowering operator

    H₀ = 0.5 * σ_z
    ℒ = liouvillian(H₀, [c_op]; convention = :TDSE)

    # Start from |0⟩⟨0|
    ρ₀ = ComplexF64[1 0; 0 0]
    ρ⃗₀ = reshape(ρ₀, :)

    tlist = collect(range(0, 10.0, length = 201))

    ρ⃗_out = propagate(
        ρ⃗₀,
        (ℒ,),
        tlist;
        method = ExponentialUtilities,
        expv_kwargs = (; ishermitian = false)
    )

    ρ_out = reshape(ρ⃗_out, 2, 2)

    # Density matrix must have trace 1 and be positive semidefinite
    @test tr(ρ_out) ≈ 1.0 atol = 1e-9
    @test all(eigvals(Hermitian(ρ_out)) .>= -1e-10)

    # Compare with exact solution via matrix exponential
    L_full = Array(ℒ)
    U = exp(-1im * L_full * tlist[end])
    ρ⃗_expected = U * ρ⃗₀
    @test norm(ρ⃗_out - ρ⃗_expected) < 1e-9
end


@testset "Liouvillian backward propagation" begin
    σ_z = ComplexF64[1 0; 0 -1]
    γ = 0.01
    c_op = √γ * ComplexF64[0 1; 0 0]

    H₀ = 0.5 * σ_z
    ℒ = liouvillian(H₀, [c_op]; convention = :TDSE)

    ρ₀ = ComplexF64[1 0; 0 0]
    ρ⃗₀ = reshape(ρ₀, :)

    tlist = collect(range(0, 1.0, length = 51))

    ρ⃗_fw = propagate(
        ρ⃗₀,
        (ℒ,),
        tlist;
        method = ExponentialUtilities,
        expv_kwargs = (; ishermitian = false)
    )

    ρ⃗_bw = propagate(
        ρ⃗_fw,
        (ℒ,),
        tlist;
        method = ExponentialUtilities,
        backward = true,
        expv_kwargs = (; ishermitian = false)
    )

    # For non-unitary dynamics, backward propagation won't exactly recover ρ₀,
    # but applying exp(-iLdt) backward should still be consistent with exp(+iLdt)
    # Compare with exact matrix exponential backward propagation
    L_full = Array(ℒ)
    U_bw = exp(+1im * L_full * tlist[end])
    ρ⃗_bw_expected = U_bw * ρ⃗_fw
    @test norm(ρ⃗_bw - ρ⃗_bw_expected) < 1e-9
end


@testset "Time-dependent generator" begin
    rng = StableRNG(20260216)
    Ψ₀ = ComplexF64[1, 0]
    H₀ = ComplexF64[0 0; 0 1]
    H₁ = ComplexF64[0 1; 1 0]

    tlist = collect(range(0, 10.0, length = 201))
    ϵ = 0.2 * ones(length(tlist))

    generator = hamiltonian(H₀, (H₁, ϵ))

    Ψ_out_expv = propagate(Ψ₀, generator, tlist; method = ExponentialUtilities)
    Ψ_out_exp = propagate(Ψ₀, generator, tlist; method = :expprop)

    @test norm(Ψ_out_expv - Ψ_out_exp) < 1e-9
end

# TODO: test a GradGen propagation
