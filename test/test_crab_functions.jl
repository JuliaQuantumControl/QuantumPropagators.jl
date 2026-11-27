# SPDX-FileCopyrightText: © 2021 Michael Goerz <mail@michaelgoerz.net>
#
# SPDX-License-Identifier: MIT

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
    @test isempty(ϵ.cos_indices)
    @test ϵ.sin_indices == 1:N
    @test (ϵ(-T) + ϵ(T)) ≈ 0.0

    ϵ = F(N; parity = :odd, guess = (t->1.0), scale_guess = true, max_frequency = 2.0, rng)
    @test isempty(ϵ.cos_indices)
    @test ϵ.sin_indices == 2:(N+1)
    if F == CRABFunction
        @test length(ϵ.parameters) == N + 1
    elseif F == VariedFrequencyCRABFunction
        @test length(ϵ.parameters) == 2N + 1
    end

    ϵ = F(N; parity = :even, max_frequency = 2.0, rng)
    @test ϵ.cos_indices == 1:N
    @test isempty(ϵ.sin_indices)
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
    @test ϵ.cos_indices == 1:N
    @test ϵ.sin_indices == (N+1):2N
    @test ϵ(-T) ≠ ϵ(T)

    ϵ = F(
        N;
        parity = :evenodd,
        guess = (t->1.0),
        scale_guess = true,
        max_frequency = 2.0,
        rng
    )
    @test ϵ.cos_indices == 2:(N+1)
    @test ϵ.sin_indices == (N+2):(2N+1)
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


# The CRAB formula, evaluated from the documented layout of the parameters
function crab_reference(ϵ, t)
    N = length(ϵ.frequencies)
    p = copy(ϵ.parameters)
    c₀ = ϵ.scale_guess ? popfirst!(p) : 1.0
    c_plus = (ϵ.parity == :odd) ? zeros(N) : splice!(p, 1:N)
    c_minus = (ϵ.parity == :even) ? zeros(N) : splice!(p, 1:N)
    r = (ϵ isa VariedFrequencyCRABFunction) ? splice!(p, 1:N) : ones(N)
    @assert isempty(p)
    ω = ϵ.frequencies
    S = isnothing(ϵ.shape) ? 1.0 : ϵ.shape(t)
    g = isnothing(ϵ.guess) ? 0.0 : ϵ.guess(t)
    return c₀ * g + S * sum(c_plus .* cos.(r .* ω .* t) .+ c_minus .* sin.(r .* ω .* t))
end

function allocations(f, t)
    f(t)  # compile
    return @allocated f(t)
end


@testset "$(nameof(F)): evaluation" for F in (CRABFunction, VariedFrequencyCRABFunction)

    rng = StableRNG(3167069968)
    N = 5
    T = 1.0
    tlist = collect(range(0, T; length = 21))
    S(t) = flattop(t; T, t_rise = 0.3T)
    guess(t) = 0.2 + t

    for parity in (:evenodd, :even, :odd), scale_guess in (true, false)
        ϵ = F(N; max_frequency = 20.0, shape = S, guess, scale_guess, parity, rng)
        ϵ.parameters .= 1.0 .- 2.0 .* rand(rng, length(ϵ.parameters))
        if F == VariedFrequencyCRABFunction
            @test ϵ.r_indices == (length(ϵ.parameters)-N+1):length(ϵ.parameters)
            ϵ.parameters[ϵ.r_indices] .= 0.5 .+ rand(rng, N)
        end
        @test ϵ.(tlist) ≈ crab_reference.(Ref(ϵ), tlist)
        @test allocations(ϵ, 0.5) == 0
    end

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
