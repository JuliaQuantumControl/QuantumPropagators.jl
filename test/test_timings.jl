# SPDX-FileCopyrightText: © 2021 Michael Goerz <mail@michaelgoerz.net>
#
# SPDX-License-Identifier: MIT

using Test
using QuantumPropagators
using QuantumPropagators: hamiltonian
using QuantumPropagators.Arnoldi: arnoldi!
using QuantumPropagators.SpectralRange: specrange
using QuantumControlTestUtils.RandomObjects:
    random_matrix, random_state_vector, random_dynamic_generator
using StableRNGs: StableRNG
using TimerOutputs
using ExponentialUtilities


@testset "Capture timings" begin

    enabled = QuantumPropagators.timings_enabled()
    @test enabled ≡ false

    QuantumPropagators.enable_timings()

    enabled = QuantumPropagators.timings_enabled()
    @test enabled ≡ true

    N = 10
    tlist = collect(range(0, 10, length = 101))
    rng = StableRNG(677918056)
    Ψ = random_state_vector(N; rng)
    Ĥ = hamiltonian(random_dynamic_generator(N, tlist; rng)...)
    propagator = init_prop(Ψ, Ĥ, tlist; method = :cheby)
    for interval = 1:(length(tlist)-1)
        prop_step!(propagator)
    end

    @test TimerOutputs.ncalls(
        propagator.wrk.timing_data["prop_step!"]["matrix-vector product"]
    ) > 200

    enabled = QuantumPropagators.timings_enabled()
    @test enabled ≡ true

    QuantumPropagators.disable_timings()
    enabled = QuantumPropagators.timings_enabled()
    @test enabled ≡ false

end


@testset "Capture timings in package extensions" begin

    # https://github.com/JuliaQuantumControl/QuantumPropagators.jl/issues/121

    extmod =
        Base.get_extension(QuantumPropagators, :QuantumPropagatorsExponentialUtilitiesExt)
    @test !isnothing(extmod)

    @test QuantumPropagators.enable_timings() ≡ true

    N = 10
    tlist = collect(range(0, 10, length = 101))
    rng = StableRNG(110814039)
    Ψ = random_state_vector(N; rng)
    Ĥ = hamiltonian(random_dynamic_generator(N, tlist; rng)...)
    propagator = init_prop(Ψ, Ĥ, tlist; method = ExponentialUtilities)
    for interval = 1:(length(tlist)-1)
        prop_step!(propagator)
    end
    @test TimerOutputs.ncalls(propagator.timing_data["prop_step!"]) == length(tlist) - 1
    @test TimerOutputs.ncalls(propagator.timing_data["prop_step!"]["expv"]) ==
          length(tlist) - 1

    # `timings_enabled` must detect a loaded extension without timings
    TimerOutputs.disable_debug_timings(extmod)
    @test QuantumPropagators.timings_enabled() ≡ false

    @test QuantumPropagators.disable_timings() ≡ false

end


@testset "Arnoldi timings without global fallback" begin

    # https://github.com/JuliaQuantumControl/QuantumPropagators.jl/issues/131

    @test !isdefined(QuantumPropagators.Arnoldi, :_DEFAULT_TIMING_DATA)

    @test QuantumPropagators.enable_timings() ≡ true

    N = 10
    m = 5
    rng = StableRNG(3261956410)
    H = random_matrix(N; hermitian = true, rng)
    Ψ = random_state_vector(N; rng)

    Hess1 = zeros(ComplexF64, m + 1, m + 1)
    q1 = [similar(Ψ) for _ = 1:(m+1)]
    timing_data = TimerOutput()
    @test arnoldi!(Hess1, q1, m, Ψ, H; _timing_data = timing_data) == m
    @test TimerOutputs.ncalls(timing_data["matrix-vector product"]) == m

    Hess2 = zeros(ComplexF64, m + 1, m + 1)
    q2 = [similar(Ψ) for _ = 1:(m+1)]
    @test arnoldi!(Hess2, q2, m, Ψ, H) == m
    @test Hess2 == Hess1

    E_min, E_max = specrange(H; method = :arnoldi)
    @test E_min < E_max

    @test QuantumPropagators.disable_timings() ≡ false

end
