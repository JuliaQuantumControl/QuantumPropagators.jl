# Test routines related to the spectral radius / range

using Test
using LinearAlgebra
using StableRNGs: StableRNG
using QuantumPropagators
using QuantumPropagators: Cheby
using QuantumPropagators.Controls: get_controls, evaluate
using QuantumPropagators.SpectralRange
using QuantumControlTestUtils.RandomObjects:
    random_matrix, random_state_vector, random_dynamic_generator


@testset "Ritz values (non-Herm)" begin

    SHOWVALS = false
    rng = StableRNG(2351042865)

    N = 1000
    ρ = 10  # spectral radius
    X = random_matrix(N; spectral_radius = ρ, rng)
    Ψ = random_state_vector(N)
    # a non_Hermitian X requires significantly more iterations to converge than
    # the default
    ritzvals = QuantumPropagators.SpectralRange.ritzvals(X, Ψ, 180, 200; prec = 1e-5)
    evals = eigvals(X)

    E_min = abs(evals[1])
    Δ_min = abs(evals[1] - ritzvals[1])
    SHOWVALS && @show evals[1], ritzvals[1]
    ϵ_min = Δ_min / E_min
    SHOWVALS && @show E_min, Δ_min, ϵ_min
    @test ϵ_min < 0.01

    E_max = abs(evals[end])
    Δ_max = abs(evals[end] - ritzvals[end])
    SHOWVALS && @show evals[end], ritzvals[end]
    ϵ_max = Δ_max / E_max
    SHOWVALS && @show E_max, Δ_max, ϵ_max
    @test ϵ_max < 0.01

    # Note: E_min and E_max should be roughly equal to ρ (Girko-Ginibri
    # circular law)

end

@testset "Ritz values (Hermitian)" begin

    SHOWVALS = false
    rng = StableRNG(2351042866)

    N = 1000
    ρ = 10  # spectral radius
    H = random_matrix(N; spectral_radius = ρ, hermitian = true, rng)
    Ψ = random_state_vector(N)
    # for a Hermitian H, lower default (cf. specrange) should be ok
    ritzvals = QuantumPropagators.SpectralRange.ritzvals(H, Ψ, 20, 60; prec = 1e-3)
    evals = eigvals(H)

    E_min = abs(evals[1])
    Δ_min = abs(evals[1] - ritzvals[1])
    SHOWVALS && @show evals[1], ritzvals[1]
    ϵ_min = Δ_min / E_min
    SHOWVALS && @show E_min, Δ_min, ϵ_min
    @test ϵ_min < 0.02

    E_max = abs(evals[end])
    Δ_max = abs(evals[end] - ritzvals[end])
    SHOWVALS && @show evals[end], ritzvals[end]
    ϵ_max = Δ_max / E_max
    SHOWVALS && @show E_max, Δ_max, ϵ_max
    @test ϵ_max < 0.02

    # Note: E_min and E_max should be roughly equal to ρ (Girko-Ginibri
    # circular law)

end


@testset "Spectral range" begin

    SHOWVALS = false
    rng = StableRNG(2351042867)

    N = 1000
    ρ = 10  # spectral radius
    density = 0.1
    H = random_matrix(N; hermitian = true, spectral_radius = ρ, density, rng)

    evals = eigvals(Array(H))
    SHOWVALS && @show evals[1], evals[end]

    prec = 1e-4
    E_min, E_max = specrange(H, method = :arnoldi, prec = prec)
    @test E_min isa Float64
    @test E_max isa Float64
    SHOWVALS && @show E_min, E_max

    Δ = evals[end] - evals[1]
    SHOWVALS && @show Δ
    δmin = E_min - evals[1]
    δmax = E_max - evals[end]
    SHOWVALS && @show δmin, δmax
    @test evals[1] - 0.05 * Δ <= E_min <= evals[1]
    @test evals[end] <= E_max < evals[end] + 0.05 * Δ

    E_min, E_max = specrange(H, method = :diag)
    @test E_min isa Float64
    @test E_max isa Float64
    SHOWVALS && @show E_min, E_max
    @test abs(evals[1] - E_min) < 1e-12
    @test abs(evals[end] - E_max) < 1e-12

    @test_throws UndefKeywordError begin
        specrange(H, method = :manual)
    end
    @test_throws UndefKeywordError begin
        specrange(H, method = :manual, E_min = -1.0)
    end
    @test_throws UndefKeywordError begin
        specrange(H, method = :manual, E_max = -1.0)
    end
    E_min, E_max = specrange(H, method = :manual, E_min = -10, E_max = 10)
    @test E_min isa Float64
    @test E_max isa Float64
    SHOWVALS && @show E_min, E_max
    @test E_min == -10.0
    @test E_max == 10.0

    E_min, E_max = specrange(H)  # method=:auto
    @test E_min isa Float64
    @test E_max isa Float64
    SHOWVALS && @show E_min, E_max
    @test evals[1] - 0.05 * Δ <= E_min <= evals[1]
    @test evals[end] <= E_max < evals[end] + 0.05 * Δ

    E_min, E_max = specrange(H; E_min = -10, E_max = 10)  # method=:auto => :manual
    @test E_min isa Float64
    @test E_max isa Float64
    SHOWVALS && @show E_min, E_max
    @test E_min == -10.0
    @test E_max == 10.0

end


@testset "Chebychev specrange" begin

    SHOWVALS = false
    rng = StableRNG(2351042868)

    N = 1000
    ρ = 10  # spectral radius
    density = 0.1
    tlist = collect(range(0, 10; length = 501))
    Ψ = random_state_vector(N)
    H_of_t = random_dynamic_generator(
        N,
        tlist;
        hermitian = true,
        spectral_envelope = ρ,
        density,
        rng
    )
    control = get_controls(H_of_t)[1]

    H_m1 = evaluate(H_of_t; vals_dict = IdDict(control => -1.0))
    H_0 = evaluate(H_of_t; vals_dict = IdDict(control => 0.0))
    H_p1 = evaluate(H_of_t; vals_dict = IdDict(control => 1.0))

    evals_m1 = eigvals(Array(H_m1))
    evals_0 = eigvals(Array(H_0))
    evals_p1 = eigvals(Array(H_p1))

    SHOWVALS && @show evals_m1[1], evals_m1[end]
    SHOWVALS && @show evals_0[1], evals_0[end]
    SHOWVALS && @show evals_p1[1], evals_p1[end]

    propagator = init_prop(Ψ, H_of_t, tlist; method = Cheby, rng)
    # we passed `rng` to make the default specrad_method=:arnoldi deterministic
    SHOWVALS && @info "Initialized with defaults" propagator.wrk.E_min propagator.wrk.Δ
    @test -10.2 < propagator.wrk.E_min < -10.1
    @test 20.1 < propagator.wrk.Δ < 20.3

    E_min = -10
    E_max = 10
    propagator = init_prop(Ψ, H_of_t, tlist; method = Cheby, E_min, E_max)
    SHOWVALS &&
        @info "Initialized with manual" E_min E_max propagator.wrk.E_min propagator.wrk.Δ
    @test propagator.wrk.E_min ≈ -10.1
    @test propagator.wrk.Δ ≈ 20.2

    propagator = init_prop(
        Ψ,
        H_of_t,
        tlist;
        method = Cheby,
        E_min,
        E_max,
        specrange_method = :manual,
        specrange_buffer = 0.1
    )
    SHOWVALS &&
        @info "Initialized with specrange_buffer=0.1" propagator.wrk.E_min propagator.wrk.Δ
    @test propagator.wrk.E_min ≈ -11.0
    @test propagator.wrk.Δ ≈ 22.0

    propagator = init_prop(
        Ψ,
        H_of_t,
        tlist;
        method = Cheby,
        specrange_method = :diag,
        specrange_buffer = 0.0,
        control_ranges = IdDict(control => (-1, 1))
    )
    SHOWVALS &&
        @info "Initialized with specrange_buffer=0" propagator.wrk.E_min propagator.wrk.Δ
    @test propagator.wrk.E_min ≈ min(evals_m1[1], evals_0[1], evals_p1[1])
    @test propagator.wrk.Δ ≈
          max(evals_m1[end], evals_0[end], evals_p1[end]) - propagator.wrk.E_min

end
