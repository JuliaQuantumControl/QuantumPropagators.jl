using Test
using LinearAlgebra
using QuantumPropagators

@testset "random Hermitian state" begin

    # Input
    N = 1000

    X = rand(ComplexF64, (N, N))
    H = Hermitian(X)

    dt = 0.5

    Ψ₀ = rand(ComplexF64, N)
    Ψ₀ ./= norm(Ψ₀)

    @test norm(Ψ₀) ≈ 1

    # Expected
    U = exp(-im * H * dt)
    Ψ_out_expected = U * Ψ₀
    @test norm(Ψ_out_expected) ≈ 1

    Ψ = copy(Ψ₀)
    wrk = ExpPropWrk(Ψ₀)
    expprop!(Ψ, H, dt, wrk)
    Ψ_out = copy(Ψ)
    @test norm(Ψ_out) ≈ 1

    # Comparison (should be exact to machine precision)
    @test norm(Ψ_out - Ψ_out_expected) ≈ 0
end

