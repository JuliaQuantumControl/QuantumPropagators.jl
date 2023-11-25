# Test routines related the the spectral radius

using Test
using LinearAlgebra
using QuantumPropagators
using QuantumPropagators.SpectralRange
using QuantumControlTestUtils.RandomObjects: random_matrix, random_state_vector


@testset "Ritz values (non-Herm)" begin

    SHOWVALS = false

    N = 1000
    ρ = 10  # spectral radius
    X = random_matrix(N; spectral_radius=ρ)
    Ψ = random_state_vector(N)
    # a non_Hermitian X requires significantly more iterations to converge than
    # the default
    ritzvals = QuantumPropagators.SpectralRange.ritzvals(X, Ψ, 180, 200; prec=1e-5)
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

    N = 1000
    ρ = 10  # spectral radius
    H = random_matrix(N; spectral_radius=ρ, hermitian=true)
    Ψ = random_state_vector(N)
    # for a Hermitian H, lower default (cf. specrange) should be ok
    ritzvals = QuantumPropagators.SpectralRange.ritzvals(H, Ψ, 20, 60; prec=1e-3)
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

    N = 1000
    ρ = 10  # spectral radius
    density = 0.1
    H = random_matrix(N; hermitian=true, spectral_radius=ρ, density)

    evals = eigvals(Array(H))
    SHOWVALS && @show evals[1], evals[end]

    prec = 1e-4
    E_min, E_max = specrange(H, method=:arnoldi, prec=prec)
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

    E_min, E_max = specrange(H, method=:diag)
    @test E_min isa Float64
    @test E_max isa Float64
    SHOWVALS && @show E_min, E_max
    @test abs(evals[1] - E_min) < 1e-12
    @test abs(evals[end] - E_max) < 1e-12

end
