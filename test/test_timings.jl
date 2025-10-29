using Test
using QuantumPropagators
using QuantumControlTestUtils.RandomObjects: random_state_vector, random_dynamic_generator
using StableRNGs: StableRNG
using TimerOutputs


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
    Ĥ = random_dynamic_generator(N, tlist; rng)
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
