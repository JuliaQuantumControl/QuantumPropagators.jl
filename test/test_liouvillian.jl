using Test
using QuantumPropagators
using QuantumPropagators.Generators
using QuantumPropagators.Controls
using LinearAlgebra
using Distributions
using QuantumControlTestUtils.RandomObjects: random_matrix, random_state_vector
using SparseArrays


@testset "TLS dissipation" begin

    function ket(i)
        Ψ = zeros(ComplexF64, 2)
        Ψ[i+1] = 1.0
        return Ψ
    end

    bra(i) = adjoint(ket(i))

    ketbra(i, j) = ket(i) * bra(j)

    γ₁ = 0.5
    γ₂ = 0.2

    Âᵧ₁ = √γ₁ * ketbra(0, 1)
    Âᵧ₂ = √(2γ₂) * ketbra(1, 1)

    Ψ₀ = (ket(0) + ket(1)) / √2
    ρ⃗₀ = reshape(Ψ₀ * Ψ₀', :)
    ℒ_op = liouvillian(nothing, [Âᵧ₁, Âᵧ₂]; convention=:TDSE)
    ℒ = Array(ℒ_op)

    T = 1.0
    tlist = [0.0, T]

    ρ⃗_out = propagate(ρ⃗₀, ℒ, tlist; method=:expprop)
    ρ_out = reshape(ρ⃗_out, 2, 2)

    #! format: off
    ρ_expected = 0.5 * ComplexF64[
            (2-exp(-γ₁*T))    (exp(-(γ₁/2 + γ₂)*T));
        (exp(-(γ₁/2 + γ₂)*T))       (exp(-γ₁*T))
    ]
    #! format: on


    @test abs(1 - tr(ρ_out)) < 1e-15 # total population
    @test abs(tr(ρ_out^2)) < 1.0  # purity
    @test norm(ρ_out - ρ_expected) < 1e-15

end


@testset "LvN" begin

    N = 100

    function ket(i)
        Ψ = zeros(ComplexF64, N)
        Ψ[i] = 1.0
        return Ψ
    end

    bra(i) = adjoint(ket(i))

    ketbra(i, j) = ket(i) * bra(j)

    Ĥ₀ = random_matrix(N; hermitian=true, spectral_radius=1)
    Ĥ₁ = random_matrix(N; hermitian=true, spectral_radius=0.1)

    ϵ(t) = 1.0
    H = (Ĥ₀, (Ĥ₁, ϵ))
    Ĥ = H[1] + H[2][1] * H[2][2](0)

    Ψ₀ = random_state_vector(N)
    ρ₀ = Ψ₀ * Ψ₀'
    ρ⃗₀ = reshape(ρ₀, :)
    𝕚 = 1im

    ℒ_nodiss = liouvillian(Ĥ; convention=:LvN)
    @test ℒ_nodiss isa SparseMatrixCSC{ComplexF64,Int64}
    @test norm(𝕚 * (Ĥ * ρ₀ - ρ₀ * Ĥ) - reshape(ℒ_nodiss * ρ⃗₀, N, N)) < 1e-15

    ℒ_nodiss = liouvillian(Ĥ; convention=:TDSE)
    @test ℒ_nodiss isa SparseMatrixCSC{ComplexF64,Int64}
    @test norm((Ĥ * ρ₀ - ρ₀ * Ĥ) - reshape(ℒ_nodiss * ρ⃗₀, N, N)) < 1e-15

    γ₁ = 0.2
    decay_to_ground = [√γ₁ * ketbra(1, i) for i ∈ 2:N]

    γ₂ = 0.1
    dephasing = [√γ₂ * ketbra(i, i) for i ∈ 1:N]

    c_ops = (decay_to_ground..., dephasing...)

    L = liouvillian(Ĥ₀, c_ops; convention=:LvN)
    @test L isa SparseMatrixCSC
    ρ̇_LvN = (
        𝕚 * (Ĥ₀ * ρ₀ - ρ₀ * Ĥ₀) +
        sum([(A * ρ₀ * A' - (A' * A * ρ₀) / 2 - (ρ₀ * A' * A) / 2) for A ∈ c_ops])
    )
    ρ̇ = reshape(L * ρ⃗₀, N, N)
    @test norm(ρ̇ - ρ̇_LvN) < 1e-15

    L = liouvillian(H, c_ops; convention=:LvN)
    t = 0.0
    ℒ = L.ops[1] + L.ops[2] * L.amplitudes[1](t)
    L0 = evaluate(L, t)
    @test norm(ℒ - Array(L0)) < 1e-12

    L = liouvillian(hamiltonian(H...), c_ops; convention=:LvN)
    L0 = evaluate(L, t)
    @test norm(ℒ - Array(L0)) < 1e-12

    ρ̇_LvN = (
        𝕚 * (Ĥ * ρ₀ - ρ₀ * Ĥ) +
        sum([(A * ρ₀ * A' - (A' * A * ρ₀) / 2 - (ρ₀ * A' * A) / 2) for A ∈ c_ops])
    )

    ρ̇ = reshape(ℒ * ρ⃗₀, N, N)

    @test norm(ρ̇ - ρ̇_LvN) < 1e-15

end
