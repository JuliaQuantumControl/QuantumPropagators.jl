using Test

using QuantumControlTestUtils.RandomObjects: random_state_vector, random_dynamic_generator
using QuantumPropagators: QuantumPropagators, init_prop
using QuantumPropagators.Interfaces: check_propagator
using StableRNGs: StableRNG
using OrdinaryDiffEq: OrdinaryDiffEq
using QuantumPropagators.Shapes: flattop


@testset "Cheby Propagator Interface" begin

    N = 10
    tlist = collect(range(0, 10, length = 101))
    rng = StableRNG(677918056)
    Ψ = random_state_vector(N; rng)
    Ĥ = random_dynamic_generator(N, tlist; rng)

    propagator = init_prop(
        Ψ,
        Ĥ,
        tlist;
        method = :cheby,
        backward = false,
        inplace = true,
        verbose = false
    )

    @test check_propagator(propagator)

    propagator = init_prop(
        Ψ,
        Ĥ,
        tlist;
        method = :cheby,
        backward = false,
        inplace = false,
        verbose = false
    )

    @test check_propagator(propagator)

    propagator = init_prop(
        Ψ,
        Ĥ,
        tlist;
        method = :cheby,
        backward = true,
        inplace = true,
        verbose = false
    )

    @test check_propagator(propagator)

end


@testset "Newton Propagator Interface" begin

    N = 10
    tlist = collect(range(0, 10, length = 101))
    rng = StableRNG(677918057)
    Ψ = random_state_vector(N; rng)
    Ĥ = random_dynamic_generator(N, tlist; rng)

    propagator = init_prop(
        Ψ,
        Ĥ,
        tlist;
        method = :newton,
        backward = false,
        inplace = true,
        verbose = false
    )

    @test check_propagator(propagator)

    @test_throws ErrorException("The Newton propagator is only implemented in-place") begin
        init_prop(
            Ψ,
            Ĥ,
            tlist;
            method = :newton,
            backward = false,
            inplace = false,
            verbose = false
        )
    end

    propagator = init_prop(
        Ψ,
        Ĥ,
        tlist;
        method = :newton,
        backward = true,
        inplace = true,
        verbose = false
    )

    @test check_propagator(propagator)

end


@testset "Exp Propagator Interface" begin

    N = 10
    rng = StableRNG(677918057)
    tlist = collect(range(0, 10, length = 101))
    Ψ = random_state_vector(N; rng)
    Ĥ = random_dynamic_generator(N, tlist; rng)

    propagator = init_prop(
        Ψ,
        Ĥ,
        tlist;
        method = :expprop,
        backward = false,
        inplace = true,
        verbose = false
    )

    @test check_propagator(propagator)

    propagator = init_prop(
        Ψ,
        Ĥ,
        tlist;
        method = :expprop,
        backward = false,
        inplace = false,
        verbose = false
    )

    @test check_propagator(propagator)

    propagator = init_prop(
        Ψ,
        Ĥ,
        tlist;
        method = :expprop,
        backward = true,
        inplace = true,
        verbose = false
    )

    @test check_propagator(propagator)

end


@testset "ODE Propagator Interface (time-continuous)" begin

    N = 10
    T = 10.0
    tlist = collect(range(0, T, length = 101))
    rng = StableRNG(677918057)
    Ψ = random_state_vector(N; rng)
    amplitudes = [t -> flattop(t; T, t_rise = 0.3 * T),]
    Ĥ = random_dynamic_generator(N, tlist; amplitudes, rng)

    propagator = init_prop(
        Ψ,
        Ĥ,
        tlist;
        method = OrdinaryDiffEq,
        backward = false,
        inplace = true,
        verbose = false
    )
    @test check_propagator(propagator)

    propagator = init_prop(
        Ψ,
        Ĥ,
        tlist;
        method = OrdinaryDiffEq,
        backward = false,
        inplace = false,
        verbose = false
    )
    @test check_propagator(propagator)

    propagator = init_prop(
        Ψ,
        Ĥ,
        tlist;
        method = OrdinaryDiffEq,
        backward = true,
        inplace = true,
        verbose = false
    )
    @test check_propagator(propagator)

end


@testset "ODE Propagator Interface (PWC)" begin

    N = 10
    T = 10.0
    tlist = collect(range(0, T, length = 101))
    rng = StableRNG(677918057)
    Ψ = random_state_vector(N; rng)
    Ĥ = random_dynamic_generator(N, tlist; rng)

    propagator = init_prop(
        Ψ,
        Ĥ,
        tlist;
        method = OrdinaryDiffEq,
        pwc = true,
        backward = false,
        inplace = true,
        verbose = false
    )
    @test check_propagator(propagator)

    propagator = init_prop(
        Ψ,
        Ĥ,
        tlist;
        method = OrdinaryDiffEq,
        pwc = true,
        backward = false,
        inplace = false,
        verbose = false
    )
    @test check_propagator(propagator)

    propagator = init_prop(
        Ψ,
        Ĥ,
        tlist;
        method = OrdinaryDiffEq,
        pwc = true,
        backward = true,
        inplace = true,
        verbose = false
    )
    @test check_propagator(propagator)

end


# Minimal struct for testing the base set_t! (all built-in propagators override it)
mutable struct _SimpleProp
    tlist::Vector{Float64}
    t::Float64
end

# Minimal AbstractPropagator that ignores piecewise/pwc, for testing enforcement errors
struct _NonPiecewiseProp <: QuantumPropagators.AbstractPropagator end

function QuantumPropagators.init_prop(state, generator, tlist, ::Val{:_non_piecewise}; _...)
    return _NonPiecewiseProp()
end


@testset "Propagator property access" begin

    N = 10
    tlist = collect(range(0, 10, length = 101))
    rng = StableRNG(677918059)
    Ψ = random_state_vector(N; rng)
    Ĥ = random_dynamic_generator(N, tlist; rng)
    propagator = init_prop(Ψ, Ĥ, tlist; method = :cheby)

    @test_throws ErrorException propagator.generator
    @test_throws ErrorException (propagator.state = Ψ)
    @test_throws ErrorException (propagator.tlist = tlist)
    @test_throws ErrorException (propagator.t = 0.0)
    @test_throws ErrorException (propagator.generator = Ĥ)

    new_params = copy(propagator.parameters)
    propagator.parameters = new_params
    @test propagator.parameters === new_params

    pub = propertynames(propagator)
    @test pub == (:state, :tlist, :t, :parameters, :backward, :inplace)

    priv = propertynames(propagator, true)
    @test length(priv) >= length(pub)
    @test all(p ∈ priv for p ∈ pub)

end


@testset "init_prop piecewise enforcement" begin

    N = 10
    tlist = collect(range(0, 10, length = 101))
    rng = StableRNG(677918060)
    Ψ = random_state_vector(N; rng)
    Ĥ = random_dynamic_generator(N, tlist; rng)

    # piecewise=true with a PWC method works (ChebyPropagator isa PiecewisePropagator)
    propagator = init_prop(Ψ, Ĥ, tlist; method = :cheby, piecewise = true)
    @test propagator isa QuantumPropagators.PiecewisePropagator

    # piecewise=true with a non-piecewise method throws
    @test_throws ErrorException init_prop(
        Ψ,
        Ĥ,
        tlist;
        method = :_non_piecewise,
        piecewise = true
    )

    # pwc=true with a non-piecewise-constant method throws
    @test_throws ErrorException init_prop(
        Ψ,
        Ĥ,
        tlist;
        method = :_non_piecewise,
        pwc = true
    )

end


@testset "init_prop unknown method" begin

    tlist = collect(range(0, 10, length = 101))
    rng = StableRNG(677918061)
    Ψ = random_state_vector(10; rng)
    Ĥ = random_dynamic_generator(10, tlist; rng)

    @test_throws ArgumentError init_prop(Ψ, Ĥ, tlist; method = 42)

end


@testset "set_t! interior and snapping" begin

    tlist = collect(range(0, 10, length = 101))
    prop = _SimpleProp(tlist, tlist[1])

    # Interior time value (the t <= tlist[1] and t >= tlist[end] branches are
    # covered by reinit_prop! calls, but the searchsortedfirst branch is not)
    QuantumPropagators.set_t!(prop, tlist[10])
    @test prop.t == tlist[10]

    # A value not on the grid triggers the snapping warning
    t_off = (tlist[10] + tlist[11]) / 2
    @test_logs (:warn, r"Snapping") QuantumPropagators.set_t!(prop, t_off)
    @test prop.t == tlist[11]

end


@testset "_get_uniform_dt non-uniform" begin

    tlist_nonuniform = [0.0, 1.0, 3.0, 6.0]

    @test QuantumPropagators._get_uniform_dt(tlist_nonuniform) === nothing

    @test_logs (:warn, r"Non-uniform") QuantumPropagators._get_uniform_dt(
        tlist_nonuniform;
        warn = true
    )

end
