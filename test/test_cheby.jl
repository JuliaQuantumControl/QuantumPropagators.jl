module TestCheby
    using Test
    using LinearAlgebra
    using QuantumPropagators

    @testset "Cheby" begin

        precision = 1e-10

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
        @test norm(U * U'  - one(U)) < precision  # U is unitary
        Ψ_out_expected = U * Ψ₀
        @test norm(Ψ_out_expected) ≈ 1

        # Cheby result
        evals = eigvals(H)

        E_min = evals[1]
        Δ = evals[end] - evals[1]

        a = cheby_coeffs(Δ, dt)
        b = zeros(20);
        cheby_coeffs!(b, Δ, dt)
        @test b[1:length(b)] ≈ a[1:length(b)]

        Ψ = copy(Ψ₀)
        wrk = ChebyWrk(Ψ₀, Δ, E_min, dt)
        cheby!(Ψ, H, dt, wrk)
        Ψ_out = copy(Ψ)

        # Comparison
        @test norm(Ψ_out - Ψ_out_expected) < precision

    end

end
