using Test

using QuantumControlTestUtils.RandomObjects: random_state_vector, random_dynamic_generator
using QuantumPropagators: init_prop
using QuantumPropagators.Interfaces: check_propagator
using StableRNGs: StableRNG
using OrdinaryDiffEq: OrdinaryDiffEq
using QuantumPropagators.Shapes: flattop


@testset "Cheby Propagator Interface" begin

    N = 10
    tlist = collect(range(0, 10, length=101))
    rng = StableRNG(677918056)
    Ψ = random_state_vector(N; rng)
    Ĥ = random_dynamic_generator(N, tlist; rng)

    propagator =
        init_prop(Ψ, Ĥ, tlist; method=:cheby, backward=false, inplace=true, verbose=false)

    @test check_propagator(propagator)

    propagator =
        init_prop(Ψ, Ĥ, tlist; method=:cheby, backward=false, inplace=false, verbose=false)

    @test check_propagator(propagator)

    propagator =
        init_prop(Ψ, Ĥ, tlist; method=:cheby, backward=true, inplace=true, verbose=false)

    @test check_propagator(propagator)

end


@testset "Newton Propagator Interface" begin

    N = 10
    tlist = collect(range(0, 10, length=101))
    rng = StableRNG(677918057)
    Ψ = random_state_vector(N; rng)
    Ĥ = random_dynamic_generator(N, tlist; rng)

    propagator =
        init_prop(Ψ, Ĥ, tlist; method=:newton, backward=false, inplace=true, verbose=false)

    @test check_propagator(propagator)

    @test_throws ErrorException("The Newton propagator is only implemented in-place") begin
        init_prop(
            Ψ,
            Ĥ,
            tlist;
            method=:newton,
            backward=false,
            inplace=false,
            verbose=false
        )
    end

    propagator =
        init_prop(Ψ, Ĥ, tlist; method=:newton, backward=true, inplace=true, verbose=false)

    @test check_propagator(propagator)

end


@testset "Exp Propagator Interface" begin

    N = 10
    rng = StableRNG(677918057)
    tlist = collect(range(0, 10, length=101))
    Ψ = random_state_vector(N; rng)
    Ĥ = random_dynamic_generator(N, tlist; rng)

    propagator = init_prop(
        Ψ,
        Ĥ,
        tlist;
        method=:expprop,
        backward=false,
        inplace=true,
        verbose=false
    )

    @test check_propagator(propagator)

    propagator = init_prop(
        Ψ,
        Ĥ,
        tlist;
        method=:expprop,
        backward=false,
        inplace=false,
        verbose=false
    )

    @test check_propagator(propagator)

    propagator =
        init_prop(Ψ, Ĥ, tlist; method=:expprop, backward=true, inplace=true, verbose=false)

    @test check_propagator(propagator)

end


@testset "ODE Propagator Interface (time-continuous)" begin

    N = 10
    T = 10.0
    tlist = collect(range(0, T, length=101))
    rng = StableRNG(677918057)
    Ψ = random_state_vector(N; rng)
    amplitudes = [t -> flattop(t; T, t_rise=0.3 * T),]
    Ĥ = random_dynamic_generator(N, tlist; amplitudes, rng)

    propagator = init_prop(
        Ψ,
        Ĥ,
        tlist;
        method=OrdinaryDiffEq,
        backward=false,
        inplace=true,
        verbose=false
    )
    @test check_propagator(propagator)

    propagator = init_prop(
        Ψ,
        Ĥ,
        tlist;
        method=OrdinaryDiffEq,
        backward=false,
        inplace=false,
        verbose=false
    )
    @test check_propagator(propagator)

    propagator = init_prop(
        Ψ,
        Ĥ,
        tlist;
        method=OrdinaryDiffEq,
        backward=true,
        inplace=true,
        verbose=false
    )
    @test check_propagator(propagator)

end


@testset "ODE Propagator Interface (PWC)" begin

    N = 10
    T = 10.0
    tlist = collect(range(0, T, length=101))
    rng = StableRNG(677918057)
    Ψ = random_state_vector(N; rng)
    Ĥ = random_dynamic_generator(N, tlist; rng)

    propagator = init_prop(
        Ψ,
        Ĥ,
        tlist;
        method=OrdinaryDiffEq,
        pwc=true,
        backward=false,
        inplace=true,
        verbose=false
    )
    @test check_propagator(propagator)

    propagator = init_prop(
        Ψ,
        Ĥ,
        tlist;
        method=OrdinaryDiffEq,
        pwc=true,
        backward=false,
        inplace=false,
        verbose=false
    )
    @test check_propagator(propagator)

    propagator = init_prop(
        Ψ,
        Ĥ,
        tlist;
        method=OrdinaryDiffEq,
        pwc=true,
        backward=true,
        inplace=true,
        verbose=false
    )
    @test check_propagator(propagator)

end
