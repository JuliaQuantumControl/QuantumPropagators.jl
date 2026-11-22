# SPDX-FileCopyrightText: © 2021 Michael Goerz <mail@michaelgoerz.net>
#
# SPDX-License-Identifier: MIT

using Test
using LinearAlgebra
using QuantumPropagators
using QuantumPropagators.Cheby
using QuantumPropagators.Amplitudes: ShapedAmplitude
using QuantumPropagators.Controls: evaluate

@testset "random state" begin

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
    @test norm(U * U' - one(U)) < precision  # U is unitary
    Ψ_out_expected = U * Ψ₀
    @test norm(Ψ_out_expected) ≈ 1

    # Cheby result
    evals = eigvals(H)

    E_min = evals[1]
    Δ = evals[end] - evals[1]

    a = cheby_coeffs(Δ, dt)
    @test length(a) ∈ [267, 268]  # numerical fluctuations
    b = zeros(20)
    n = cheby_coeffs!(b, Δ, dt)
    @test b[1:n] ≈ a[1:n]

    Ψ = copy(Ψ₀)
    wrk = ChebyWrk(Ψ₀, Δ, E_min, dt)
    cheby!(Ψ, H, dt, wrk)
    Ψ_out = copy(Ψ)

    # Comparison
    @test norm(Ψ_out - Ψ_out_expected) < precision

end


@testset "check_normalization for non-normal generator" begin

    # A non-normal operator with the same (real) spectrum {0, 1} as Ĥ,
    # structured like a gradient generator [Ĥ μ̂; 0 Ĥ]. Its Rayleigh quotients
    # exceed the spectral range, e.g., for v = (|e⟩, |e⟩)/√2 with the top
    # eigenvector |e⟩ of Ĥ: ⟨v|G|v⟩ = E_max + ⟨e|μ̂|e⟩/2 = 1.5
    Ĥ = ComplexF64[0 0; 0 1]
    μ̂ = ComplexF64[0 0; 0 1]
    G = [Ĥ μ̂; zero(Ĥ) Ĥ]
    @test !ishermitian(G)
    @test sort(real.(eigvals(G))) ≈ [0, 0, 1, 1]

    Ψ₀ = ComplexF64[0, 1, 0, 1] / √2
    dt = 0.05
    tlist = collect(range(0, 1.0, step = dt))
    Ψ_expected = exp(-im * G * tlist[end]) * Ψ₀

    # The low-level check is only valid for Hermitian operators
    wrk = ChebyWrk(Ψ₀, 1.0, 0.0, dt)
    @test wrk.n_coeffs > 2
    @test_throws AssertionError cheby!(copy(Ψ₀), G, dt, wrk; check_normalization = true)

    # … so the propagator ignores it for a non-Hermitian generator
    kwargs = (; specrange_method = :manual, E_min = 0.0, E_max = 1.0)
    for inplace in (true, false)
        Ψ = @test_logs (:warn, r"Ignoring `check_normalization=true`") propagate(
            Ψ₀,
            G,
            tlist;
            method = Cheby,
            check_normalization = true,
            inplace,
            kwargs...
        )
        @test norm(Ψ - Ψ_expected) < 1e-10
    end

    # For a Hermitian generator, the check stays active, without a warning
    propagator = @test_logs init_prop(
        Ψ₀[3:4],
        Ĥ,
        tlist;
        method = Cheby,
        check_normalization = true,
        kwargs...
    )
    @test propagator.check_normalization
    Ψ = @test_logs propagate(
        Ψ₀[3:4],
        Ĥ,
        tlist;
        method = Cheby,
        check_normalization = true,
        kwargs...
    )
    @test norm(Ψ - exp(-im * Ĥ * tlist[end]) * Ψ₀[3:4]) < 1e-10

end


@testset "check_normalization for Hermitian sum of non-Hermitian terms" begin

    # A drive Ω â + Ω* â† in the rotating frame: `ishermitian` is `false` for
    # the `Operator`, since its terms are not Hermitian, but the sum is
    N = 6
    â = diagm(1 => sqrt.(1.0:(N-1)))
    Ĥ₀ = diagm(0 => collect(0.0:(N-1)))
    tlist = collect(range(0, 10, length = 201))
    ϵ(t) = 0.2 * sin(π * t / 10)
    Ω = ShapedAmplitude(ϵ; shape = t -> exp(0.3im))
    Ω̄ = ShapedAmplitude(ϵ; shape = t -> exp(-0.3im))
    H = hamiltonian(Ĥ₀, (â, Ω), (â', Ω̄))
    G = evaluate(H, tlist, 100)
    @test !ishermitian(G)
    @test ishermitian(Array(G))
    Ψ₀ = zeros(ComplexF64, N)
    Ψ₀[1] = 1
    propagator =
        @test_logs init_prop(Ψ₀, H, tlist; method = Cheby, check_normalization = true)
    @test propagator.check_normalization

end


@testset "check_normalization for operator without ishermitian" begin

    # Without `ishermitian`, the check stays active
    struct NoIsHermitianOperator
        H::Matrix{ComplexF64}
    end
    LinearAlgebra.mul!(y, O::NoIsHermitianOperator, x) = mul!(y, O.H, x)
    QuantumPropagators.Controls.get_controls(::NoIsHermitianOperator) = ()
    QuantumPropagators.Controls.evaluate(O::NoIsHermitianOperator, args...; kwargs...) = O
    H = NoIsHermitianOperator(ComplexF64[0 1; 1 0])
    @test_throws MethodError ishermitian(H)
    tlist = collect(range(0, 1.0, length = 11))
    propagator = @test_logs init_prop(
        ComplexF64[1, 0],
        H,
        tlist;
        method = Cheby,
        check_normalization = true,
        specrange_method = :manual,
        E_min = -1.0,
        E_max = 1.0
    )
    @test propagator.check_normalization

end
