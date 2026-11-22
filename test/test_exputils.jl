# SPDX-FileCopyrightText: © 2021 Michael Goerz <mail@michaelgoerz.net>
#
# SPDX-License-Identifier: MIT

using Test
using LinearAlgebra
using SparseArrays
using QuantumPropagators
using QuantumPropagators.Generators: liouvillian
using QuantumPropagators: Cheby
using ExponentialUtilities
using QuantumGradientGenerators: GradGenerator, GradVector
using StableRNGs: StableRNG
using QuantumControlTestUtils.RandomObjects: random_state_vector, random_matrix


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

@testset "GradGen propagation" begin
    Ψ₀ = ComplexF64[1, 0]
    H₀ = ComplexF64[0 0; 0 1]
    H₁ = ComplexF64[0 1; 1 0]

    tlist = collect(range(0, 10.0, length = 201))
    ϵ = 0.2 * ones(length(tlist))

    generator = hamiltonian(H₀, (H₁, ϵ))
    G = GradGenerator(generator)
    Ψ = GradVector(Ψ₀, 1)

    Ψ_out = propagate(
        Ψ,
        G,
        tlist;
        method = ExponentialUtilities,
        inplace = false,
        expv_kwargs = (; ishermitian = false, tol = 1e-12)
    )
    Ψ_cheby = propagate(Ψ₀, generator, tlist; method = Cheby)

    @test norm(Ψ_out.state - Ψ_cheby) < 1e-7
    @test norm(Ψ_out.grad_states[1]) > 0

    # `ExponentialUtilities.expv!` can only write to an `AbstractVector`, so
    # the in-place propagation of a `GradVector` goes through a flat buffer
    propagator = init_prop(
        GradVector(Ψ₀, 1),
        G,
        tlist;
        method = ExponentialUtilities,
        inplace = true,
        expv_kwargs = (; ishermitian = false, tol = 1e-12)
    )
    @test propagator.buffer isa Vector{ComplexF64}
    @test length(propagator.buffer) == 4
    Ψ_out_inplace = propagate(
        GradVector(Ψ₀, 1),
        G,
        tlist;
        method = ExponentialUtilities,
        inplace = true,
        expv_kwargs = (; ishermitian = false, tol = 1e-12)
    )
    @test Ψ_out_inplace isa GradVector
    @test norm(Ψ_out_inplace - Ψ_out) < 1e-12

    # Backward propagation (dt < 0) through the buffer matches out-of-place
    Ψ_bw = Dict(
        inplace => propagate(
            GradVector(Ψ_cheby, 1),
            G,
            tlist;
            method = ExponentialUtilities,
            inplace,
            backward = true,
            expv_kwargs = (; ishermitian = false, tol = 1e-12)
        ) for inplace in (true, false)
    )
    @test Ψ_bw[true] isa GradVector
    @test norm(Ψ_bw[true] - Ψ_bw[false]) < 1e-12
    @test norm(Ψ_bw[true].state - Ψ₀) < 1e-6
end


@testset "In-place propagation of state without vector interface" begin
    struct NoVectorInterfaceState
        data::Vector{ComplexF64}
    end
    tlist = collect(range(0, 1.0, length = 11))
    @test_throws ArgumentError init_prop(
        NoVectorInterfaceState(ComplexF64[1, 0]),
        ComplexF64[0 1; 1 0],
        tlist;
        method = ExponentialUtilities,
        inplace = true
    )
end


@testset "Error-estimate mode with non-Hermitian generator" begin
    rng = StableRNG(3864141668)
    N = 10
    tlist = collect(range(0, 1.0, length = 11))
    Ψ₀ = random_state_vector(N; rng)
    H = random_matrix(N; rng, hermitian = false)
    for inplace in (true, false)
        @test_throws ArgumentError(
            "`mode=:error_estimate` in `expv_kwargs` requires a Hermitian generator. Use the default `mode=:happy_breakdown` for non-Hermitian generators."
        ) init_prop(
            Ψ₀,
            H,
            tlist;
            method = ExponentialUtilities,
            inplace,
            expv_kwargs = (; mode = :error_estimate)
        )
    end
end
