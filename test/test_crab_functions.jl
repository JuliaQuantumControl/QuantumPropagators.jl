using Test

using IOCapture: IOCapture
using LinearAlgebra: norm

using StableRNGs: StableRNG
using QuantumPropagators.ParameterizedFunctions: CRABFunction, VariedFrequencyCRABFunction
using QuantumPropagators.Interfaces: check_parameterized_function
using QuantumPropagators.Shapes: flattop
using QuantumPropagators.Controls: discretize


@testset "$(nameof(F)) interface" for F in (CRABFunction, VariedFrequencyCRABFunction)

    rng = StableRNG(530280291)
    ϵ = F(10; max_frequency = 100.0, rng)
    T = 1.0
    tlist = collect(range(0, T; length = 101))
    @test check_parameterized_function(ϵ; tlist)

    rng = StableRNG(530280291)
    S(t) = flattop(t; T, t_rise = 0.1T)
    ϵ = F(10; shape = S, max_frequency = 100.0, rng)
    @test check_parameterized_function(ϵ; tlist)

    rng = StableRNG(530280291)
    ϵ = F(10; guess = (t -> 5.0), shape = S, max_frequency = 100.0, rng)
    @test check_parameterized_function(ϵ; tlist)

    rng = StableRNG(530280291)
    ϵ = F(10; guess = (t -> 5.0), shape = S, max_frequency = 100.0, parity = :odd, rng)
    @test check_parameterized_function(ϵ; tlist)

end


@testset "$(nameof(F)): no frequencies" for F in (CRABFunction, VariedFrequencyCRABFunction)

    rng = StableRNG(530280292)

    ϵ = F(10; frequencies = rand(rng, 10), rng)
    @test length(ϵ.frequencies) == 10

    c = IOCapture.capture(rethrow = Union{}, passthrough = false) do
        F(10; rng)
    end
    @test c.error
    @test contains(c.value.msg, "Did you forget to pass `max_frequency`?")

end


@testset "$(nameof(F)): parity" for F in (CRABFunction, VariedFrequencyCRABFunction)
    rng = StableRNG(530280293)

    T = -10.0
    tlist = collect(range(-T, T; length = 201))

    N = 10

    ϵ = F(N; parity = :odd, max_frequency = 2.0, rng)
    @test ϵ.i_sin == 0
    @test (ϵ(-T) + ϵ(T)) ≈ 0.0

    ϵ = F(N; parity = :odd, guess = (t->1.0), scale_guess = true, max_frequency = 2.0, rng)
    @test ϵ.i_sin == 1
    if F == CRABFunction
        @test length(ϵ.parameters) == N + 1
    elseif F == VariedFrequencyCRABFunction
        @test length(ϵ.parameters) == 2N + 1
    end

    ϵ = F(N; parity = :even, max_frequency = 2.0, rng)
    @test ϵ.i_cos == 0
    if F == CRABFunction
        @test length(ϵ.parameters) == N
    elseif F == VariedFrequencyCRABFunction
        @test length(ϵ.parameters) == 2N
    end
    @test (ϵ(-T) - ϵ(T)) ≈ 0.0

    ϵ = F(N; parity = :even, guess = (t->1.0), scale_guess = true, max_frequency = 2.0, rng)
    if F == CRABFunction
        @test length(ϵ.parameters) == N + 1
    elseif F == VariedFrequencyCRABFunction
        @test length(ϵ.parameters) == 2N + 1
    end

    ϵ = F(N; parity = :evenodd, max_frequency = 2.0, rng)
    @test ϵ.i_cos == 0
    @test ϵ.i_sin == N
    @test ϵ(-T) ≠ ϵ(T)

    ϵ = F(
        N;
        parity = :evenodd,
        guess = (t->1.0),
        scale_guess = true,
        max_frequency = 2.0,
        rng
    )
    @test ϵ.i_cos == 1
    @test ϵ.i_sin == N+1
    if F == CRABFunction
        @test length(ϵ.parameters) == 2N + 1
    elseif F == VariedFrequencyCRABFunction
        @test length(ϵ.parameters) == 3N + 1
    end
    @test ϵ(-T) ≠ ϵ(T)

    ϵ = F(N; parity = :oddeven, max_frequency = 2.0, rng)
    @test ϵ(-T) ≠ ϵ(T)

    c = IOCapture.capture(rethrow = Union{}, passthrough = false) do
        ϵ = F(N; parity = :oddevn, max_frequency = 2.0, rng)  # deliberate typo
    end
    @test c.error
    msg = "parity must be one of `:evenodd`, `:odd`, `:even`, not `:oddevn`"
    @test c.value.msg == msg

end


@testset "$(nameof(F)): rng reproducibility" for F in
                                                 (CRABFunction, VariedFrequencyCRABFunction)

    # checks that we're passing along the RNG everywhere

    T = 1.0
    tlist = collect(range(0, T; length = 101))

    rng = StableRNG(530280295)
    ϵ1 = F(10; max_frequency = 100.0, rng)

    rng = StableRNG(530280295)
    ϵ2 = F(10; max_frequency = 100.0, rng)

    @test norm(ϵ1.(tlist) .- ϵ2.(tlist)) ≈ 0.0

    rng = StableRNG(530280296)
    ϵ2 = F(10; max_frequency = 100.0, rng)
    @test norm(ϵ1.(tlist) .- ϵ2.(tlist)) > 0.0

end

@testset "$(nameof(F)): invalid frequencies" for F in
                                                 (CRABFunction, VariedFrequencyCRABFunction)

    rng = StableRNG(530280294)

    # Test when frequency length doesn't match N
    c = IOCapture.capture(rethrow = Union{}, passthrough = false) do
        F(10; frequencies = [1.0, 2.0, 3.0], rng)  # only 3 frequencies for N=10
    end
    @test c.error
    @test contains(c.value.msg, "Length of frequencies 3 must match N=10")

end


@testset "$(nameof(F)): unexpected number of parameters" for F in (
    CRABFunction,
    VariedFrequencyCRABFunction
)

    rng = StableRNG(530280295)
    c = IOCapture.capture(rethrow = Union{}, passthrough = false) do
        F(5; parameters = [1.0, 2.0], max_frequency = 10.0, rng)
    end
    @test c.error
    @test contains(c.value.msg, "Number of parameters must be")

end

@testset "Recursive guess in CRAB control" begin

    T = 1.0
    tlist = collect(range(0, T; length = 101))
    S(t) = flattop(t; T, t_rise = 0.1T)

    rng = StableRNG(530280295)
    ϵ1 = CRABFunction(10; shape = S, max_frequency = 100.0, random_amplitude = true, rng)
    ϵ2 = CRABFunction(
        10;
        shape = S,
        guess = t->1.0,
        frequencies = ϵ1.frequencies,
        parameters = ϵ1.parameters,
        scale_guess = false
    )
    @test norm((1.0 .+ ϵ1.(tlist)) .- ϵ2.(tlist)) ≈ 0.0

    rng = StableRNG(530280295)
    ϵ1 = VariedFrequencyCRABFunction(
        10;
        shape = S,
        max_frequency = 100.0,
        random_amplitude = true,
        rng
    )
    ϵ2 = VariedFrequencyCRABFunction(
        10;
        shape = S,
        guess = t->1.0,
        frequencies = ϵ1.frequencies,
        parameters = ϵ1.parameters,
        scale_guess = false
    )
    @test norm((1.0 .+ ϵ1.(tlist)) .- ϵ2.(tlist)) ≈ 0.0

    rng = StableRNG(530280295)
    ϵ1 = VariedFrequencyCRABFunction(
        10;
        shape = S,
        max_frequency = 100.0,
        random_amplitude = true,
        rng
    )
    ϵ2 = CRABFunction(
        10;
        shape = S,
        guess = ϵ1,
        random_amplitude = false,
        frequencies = ϵ1.frequencies
    )
    @test norm(ϵ1.(tlist) .- ϵ2.(tlist)) ≈ 0.0

    @test ϵ2.scale_guess
    ϵ2.parameters[1] = 0.5  # parameter for scaling guess
    @test norm(ϵ1.(tlist) .- 2 .* ϵ2.(tlist)) ≈ 0.0

    ϵ3 = CRABFunction(20; shape = S, guess = ϵ2, max_frequency = 100.0)
    @test length(ϵ3.parameters) == 41
    @test ϵ3.parameters[1] == 1.0
    @test norm(ϵ3.parameters) > 1.1

    ϵ3 = CRABFunction(20; shape = S, guess = ϵ2, max_frequency = 100.0, scale_guess = false)
    @test length(ϵ3.parameters) == 40

end


@testset "Tuneable CRAB frequency initialization" begin

    N = 10
    T = 1.0
    S(t) = flattop(t; T, t_rise = 0.1T)

    rng = StableRNG(530280298)
    frequencies = sort(100 * rand(rng, N))

    rng = StableRNG(530280299)
    ϵ1 = CRABFunction(N; shape = S, random_amplitude = true, frequencies, rng)

    rng = StableRNG(530280299)
    ϵ2 =
        VariedFrequencyCRABFunction(N; shape = S, random_amplitude = true, frequencies, rng)

    # parameters for tuneable frequencies should be initialized to 1, so the
    # pulses with variable and static frequencies should be identical.

    tlist = collect(range(0, T; length = 101))
    @test norm(ϵ1.(tlist) .- ϵ2.(tlist)) < 1e-12

end


@testset "Discretize $(nameof(F)) pulses" for F in
                                              (CRABFunction, VariedFrequencyCRABFunction)

    rng = StableRNG(530280300)
    N = 10
    T = 1.0
    S(t) = flattop(t; T, t_rise = 0.1T)
    tlist = collect(range(0, T; length = 101))

    ϵ = F(10; shape = S, max_frequency = 100.0, rng)
    pulse = discretize(ϵ, tlist)
    @test pulse isa Vector{Float64}
    @test length(pulse) == length(tlist)

    c = IOCapture.capture(rethrow = Union{}, passthrough = false) do
        F(10; shape = S, max_frequency = 100.0, rng, guess = pulse)
    end
    @test c.error
    msg = "cannot be instantiated with a vector of pulse values as a guess"
    @test contains(c.value.msg, msg)

end
