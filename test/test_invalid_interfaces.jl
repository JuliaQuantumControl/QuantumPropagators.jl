using Test
using Logging: with_logger
using QuantumControlTestUtils: QuantumTestLogger
using QuantumControlTestUtils.RandomObjects: random_dynamic_generator, random_state_vector
using StableRNGs: StableRNG
using QuantumPropagators: init_prop
using LinearAlgebra
using QuantumPropagators.Interfaces:
    check_amplitude,
    check_control,
    check_operator,
    check_generator,
    check_state,
    check_propagator


@testset "Invalid amplitude" begin

    struct InvalidAmplitude end

    tlist = collect(range(0, 10, length=101))

    ampl = InvalidAmplitude()
    test_logger = QuantumTestLogger()
    with_logger(test_logger) do
        @test check_amplitude(ampl; tlist) ≡ false
    end

    @test "get_controls(ampl)` must be defined" ∈ test_logger
    @test "all controls in `ampl` must pass `check_control`" ∈ test_logger
    @test "`substitute(ampl, replacements)` must be defined" ∈ test_logger
    @test "`evaluate(ampl, tlist, 1)` must return a Number" ∈ test_logger
    @test "`evaluate(ampl, tlist, n; vals_dict)` must be defined" ∈ test_logger

end


@testset "Invalid control" begin

    struct InvalidControl end

    tlist = collect(range(0, 10, length=101))

    control = InvalidControl()
    test_logger = QuantumTestLogger()
    with_logger(test_logger) do
        @test check_control(control; tlist) ≡ false
    end

    @test "`discretize(control, tlist)` must be defined" ∈ test_logger
    @test "`discretize_on_midpoints(control, tlist)` must be defined" ∈ test_logger

    tlist = [0.0, 1.0, 2.0, 3.0]
    control = [0.0, NaN, Inf, 0.0]
    test_logger = QuantumTestLogger()
    with_logger(test_logger) do
        @test check_control(control; tlist) ≡ false
    end
    @test "all values in `discretize(control, tlist)` must be finite" ∈ test_logger
    @test "all values in `discretize_on_midpoints(control, tlist)` must be finite" ∈
          test_logger

end


@testset "Invalid operator" begin

    struct InvalidOperator end

    state = ComplexF64[1, 0, 0, 0]
    tlist = collect(range(0, 10, length=101))

    operator = InvalidOperator()
    test_logger = QuantumTestLogger()
    with_logger(test_logger) do
        @test check_operator(operator; state, tlist) ≡ false
    end

    @test "op must not contain any controls" ∈ test_logger
    @test "`op * state` must be defined" ∈ test_logger
    @test "The 3-argument `mul!` must apply `op` to the given `state`" ∈ test_logger
    @test "The 5-argument `mul!` must apply `op` to the given `state`" ∈ test_logger
    @test "`dot(state, op, state)` must return return a number" ∈ test_logger

end


@testset "Invalid generator" begin

    struct InvalidGenerator end

    state = ComplexF64[1, 0, 0, 0]
    tlist = collect(range(0, 10, length=101))

    generator = InvalidGenerator()
    test_logger = QuantumTestLogger()
    with_logger(test_logger) do
        @test check_generator(generator; state, tlist) ≡ false
    end

    @test "`get_controls(generator)` must be defined" ∈ test_logger
    @test "`evaluate(generator, tlist, n)` must return an operator that passes `check_operator`" ∈
          test_logger
    @test "`substitute(generator, replacements)` must be defined" ∈ test_logger

end


@testset "Invalid state" begin

    struct InvalidState end
    state = InvalidState()
    test_logger = QuantumTestLogger()
    with_logger(test_logger) do
        @test check_state(state; normalized=true) ≡ false
    end
    @test "`similar(state)` must be defined" ∈ test_logger
    @test "the inner product of two states must be a complex number" ∈ test_logger
    @test "the norm of a state must be defined via the inner product" ∈ test_logger
    @test "`state + state` and `state - state` must be defined" ∈ test_logger
    @test "copy(state) must be defined" ∈ test_logger
    @test "`c * state` for a scalar `c` must be defined" ∈ test_logger
    @test "`0.0 * state` must produce a state with norm 0" ∈ test_logger
    @test "`norm(state)` must be 1" ∈ test_logger


    struct InvalidState2 end
    state = InvalidState2()
    Base.similar(::InvalidState2) = InvalidState2()
    test_logger = QuantumTestLogger()
    with_logger(test_logger) do
        @test check_state(state; normalized=true) ≡ false
    end
    @test "`copyto!(other, state)` must be defined" ∈ test_logger

    struct InvalidState3 end
    state = InvalidState3()
    Base.similar(::InvalidState3) = InvalidState3()
    Base.copyto!(a::InvalidState3, b::InvalidState3) = a
    test_logger = QuantumTestLogger()
    with_logger(test_logger) do
        @test check_state(state; normalized=true) ≡ false
    end
    @test "`similar(state)` must return a valid state" ∈ test_logger
    @test "On `similar(state)`: " ∈ test_logger


    struct InvalidState4
        Ψ
    end
    Base.similar(state::InvalidState4) = InvalidState4(similar(state.Ψ))
    Base.copy(state::InvalidState4) = InvalidState4(copy(state.Ψ))
    Base.copyto!(a::InvalidState4, b::InvalidState4) = copyto!(a.Ψ, b.Ψ)
    LinearAlgebra.dot(a::InvalidState4, b::InvalidState4) = a.Ψ ⋅ b.Ψ
    LinearAlgebra.norm(a::InvalidState4) = norm(a.Ψ)
    Base.:+(a::InvalidState4, b::InvalidState4) = InvalidState4(a.Ψ + b.Ψ)
    Base.:-(a::InvalidState4, b::InvalidState4) = InvalidState4(a.Ψ - b.Ψ)
    Base.:*(α::Number, state::InvalidState4) = InvalidState4(α * state.Ψ)
    state = InvalidState4([1im, 0])
    test_logger = QuantumTestLogger()
    with_logger(test_logger) do
        @test check_state(state; normalized=true) ≡ false
    end
    @test "`lmul!(c, state)` for a scalar `c` must be defined" ∈ test_logger
    @test "`lmul!(0.0, state)` must produce a state with norm 0" ∈ test_logger
    @test "`axpy!(c, state, other)` must be defined" ∈ test_logger

    state = [1, 0, 0, 0]
    test_logger = QuantumTestLogger()
    with_logger(test_logger) do
        @test check_state(state; normalized=true) ≡ false
    end
    @test "`state ⋅ state` must return a Complex number type" ∈ test_logger

end


@testset "Invalid propagator" begin

    struct InvalidPropagatorEmpty end

    include("invalid_propagator.jl")

    propagator = InvalidPropagatorEmpty()

    test_logger = QuantumTestLogger()
    with_logger(test_logger) do
        @test check_propagator(propagator) ≡ false
    end

    @test "does not have the required properties" ∈ test_logger

    N = 10
    tlist = collect(range(0, 100, length=1001))
    rng = StableRNG(93655235)
    Ĥ = random_dynamic_generator(N, tlist; rng)
    Ψ = random_state_vector(N; rng)

    propagator = init_prop(Ψ, Ĥ, tlist; method=:invalid_propagator_no_methods)

    test_logger = QuantumTestLogger()
    with_logger(test_logger) do
        @test check_propagator(propagator) ≡ false
    end

    @test "`prop_step!(propagator)` must be defined" ∈ test_logger
    @test "Failed to run `prop_step!(propagator)`" ∈ test_logger
    @test "`reinit_prop!` must be defined" ∈ test_logger


    propagator = init_prop(
        Ψ,
        Ĥ,
        tlist;
        method=:invalid_propagator_empty_methods,
        inplace=true,
        backward=false
    )
    test_logger = QuantumTestLogger()
    with_logger(test_logger) do
        @test check_propagator(propagator) ≡ false
    end
    @test "propagator.t ≠ propagator.tlist[begin]" ∈ test_logger
    @test "`set_t!(propagator, t)` must set propagator.t" ∈ test_logger
    @test "`prop_step!(propagator)` at final t=0.1 must return `nothing`" ∈ test_logger
    @test "`propagator.parameters` must be a dict" ∈ test_logger
    @test "`reinit_prop!(propagator, state)` must reset `propagator.t`" ∈ test_logger

    propagator = init_prop(
        Ψ,
        Ĥ,
        tlist;
        method=:invalid_propagator_empty_methods,
        inplace=false,
        backward=true
    )
    test_logger = QuantumTestLogger()
    with_logger(test_logger) do
        @test check_propagator(propagator) ≡ false
    end
    @test "propagator.t ≠ propagator.tlist[end]" ∈ test_logger
    @test "For a not-in-place propagator, the state returned by `prop_step!` must be a new object" ∈
          test_logger
    @test "`prop_step!` must advance `propagator.t` forward or backward one step on the time grid" ∈
          test_logger

    propagator = init_prop(
        Ψ,
        Ĥ,
        tlist;
        method=:invalid_random_propagator,
        inplace=true,
        backward=false
    )
    test_logger = QuantumTestLogger()
    with_logger(test_logger) do
        @test check_propagator(propagator) ≡ false
    end
    @test "For an in-place propagator, the state returned by `prop_step!` must be the `propagator.state` object" ∈
          test_logger
    @test "`set_state!(propagator, state)` for an in-place propagator must overwrite `propagator.state` in-place." ∈
          test_logger
    @test "`reinit_prop!(propagator, state)` must be idempotent" ∈ test_logger

    propagator = init_prop(
        Ψ,
        Ĥ,
        tlist;
        method=:invalid_random_propagator,
        inplace=false,
        backward=true
    )
    test_logger = QuantumTestLogger()
    with_logger(test_logger) do
        @test check_propagator(propagator) ≡ false
    end
    @test "For a not-in-place propagator, the state returned by `prop_step!` must be a new object" ∈
          test_logger
    @test "`reinit_prop!` must be defined and re-initialize the propagator" ∈ test_logger

    propagator = init_prop(
        Ψ,
        Ĥ,
        tlist;
        method=:invalid_propagator_no_state,
        inplace=true,
        backward=false
    )
    test_logger = QuantumTestLogger()
    with_logger(test_logger) do
        @test check_propagator(propagator) ≡ false
    end
    @test "prop_step! must return a valid state until time grid is exhausted" ∈ test_logger


end
