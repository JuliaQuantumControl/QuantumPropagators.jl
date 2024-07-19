using Test
using Logging: with_logger
using IOCapture: IOCapture
using QuantumControlTestUtils.RandomObjects: random_dynamic_generator, random_state_vector
using StableRNGs: StableRNG
using QuantumPropagators: QuantumPropagators, init_prop
using LinearAlgebra
import QuantumPropagators.Interfaces: supports_inplace
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
    captured = IOCapture.capture() do
        check_amplitude(ampl; tlist)
    end
    @test captured.value ≡ false

    @test contains(captured.output, "get_controls(ampl)` must be defined")
    @test contains(captured.output, "all controls in `ampl` must pass `check_control`")
    @test contains(captured.output, "`substitute(ampl, replacements)` must be defined")
    @test contains(captured.output, "`evaluate(ampl, tlist, 1)` must return a Number")
    @test contains(captured.output, "`evaluate(ampl, tlist, n; vals_dict)` must be defined")

end


@testset "Invalid control" begin

    struct InvalidControl end

    tlist = collect(range(0, 10, length=101))

    control = InvalidControl()
    captured = IOCapture.capture() do
        check_control(control; tlist)
    end
    @test captured.value ≡ false
    @test contains(captured.output, "`discretize(control, tlist)` must be defined")
    @test contains(
        captured.output,
        "`discretize_on_midpoints(control, tlist)` must be defined"
    )

    tlist = [0.0, 1.0, 2.0, 3.0]
    control = [0.0, NaN, Inf, 0.0]
    captured = IOCapture.capture() do
        check_control(control; tlist)
    end
    @test captured.value ≡ false
    @test contains(
        captured.output,
        "all values in `discretize(control, tlist)` must be finite"
    )
    @test contains(
        captured.output,
        "all values in `discretize_on_midpoints(control, tlist)` must be finite"
    )

end


@testset "Invalid operator" begin

    struct InvalidOperator end

    state = ComplexF64[1, 0, 0, 0]
    tlist = collect(range(0, 10, length=101))

    operator = InvalidOperator()
    captured = IOCapture.capture() do
        check_operator(operator; state, tlist)
    end
    @test captured.value ≡ false
    @test contains(
        captured.output,
        "The `QuantumPropagators.Interfaces.supports_inplace` method must be defined for `op`"
    )
    @test contains(captured.output, "op must not contain any controls")
    @test contains(captured.output, "`op * state` must be defined")
    @test contains(
        captured.output,
        "The 3-argument `mul!` must apply `op` to the given `state`"
    )
    @test contains(
        captured.output,
        "The 5-argument `mul!` must apply `op` to the given `state`"
    )
    @test contains(captured.output, "`dot(state, op, state)` must return return a number")

end


@testset "Invalid generator" begin

    struct InvalidGenerator end

    state = ComplexF64[1, 0, 0, 0]
    tlist = collect(range(0, 10, length=101))

    generator = InvalidGenerator()
    captured = IOCapture.capture() do
        check_generator(generator; state, tlist)
    end
    @test captured.value ≡ false

    @test contains(captured.output, "`get_controls(generator)` must be defined")
    @test contains(
        captured.output,
        "`evaluate(generator, tlist, n)` must return an operator that passes `check_operator`"
    )
    @test contains(captured.output, "`substitute(generator, replacements)` must be defined")

end


@testset "Invalid state" begin

    struct InvalidStaticState end
    state = InvalidStaticState()
    captured = IOCapture.capture() do
        check_state(state)
    end
    @test captured.value ≡ false
    @test contains(
        captured.output,
        "The `QuantumPropagators.Interfaces.supports_inplace` method must be defined for `state`"
    )

    struct InvalidState end
    QuantumPropagators.Interfaces.supports_inplace(::InvalidState) = true
    state = InvalidState()
    captured = IOCapture.capture() do
        check_state(state; normalized=true)
    end
    @test captured.value ≡ false
    @test contains(captured.output, "`similar(state)` must be defined")
    @test contains(
        captured.output,
        "the inner product of two states must be a complex number"
    )
    @test contains(
        captured.output,
        "the norm of a state must be defined via the inner product"
    )
    @test contains(captured.output, "`state + state` and `state - state` must be defined")
    @test contains(captured.output, "copy(state) must be defined")
    @test contains(captured.output, "`c * state` for a scalar `c` must be defined")
    @test contains(captured.output, "`0.0 * state` must produce a state with norm 0")
    @test contains(captured.output, "`norm(state)` must be 1")

    struct InvalidState2 end
    QuantumPropagators.Interfaces.supports_inplace(::InvalidState2) = true
    state = InvalidState2()
    Base.similar(::InvalidState2) = InvalidState2()
    captured = IOCapture.capture() do
        check_state(state; normalized=true)
    end
    @test captured.value ≡ false
    @test contains(captured.output, "`copyto!(other, state)` must be defined")

    struct InvalidState3 end
    QuantumPropagators.Interfaces.supports_inplace(::InvalidState3) = true
    state = InvalidState3()
    Base.similar(::InvalidState3) = InvalidState3()
    Base.copyto!(a::InvalidState3, b::InvalidState3) = a
    captured = IOCapture.capture() do
        check_state(state; normalized=true)
    end
    @test captured.value ≡ false
    @test contains(captured.output, "`similar(state)` must return a valid state")
    @test contains(captured.output, "On `similar(state)`: ")

    struct InvalidState4
        Ψ
    end
    QuantumPropagators.Interfaces.supports_inplace(::InvalidState4) = true
    Base.similar(state::InvalidState4) = InvalidState4(similar(state.Ψ))
    Base.copy(state::InvalidState4) = InvalidState4(copy(state.Ψ))
    Base.copyto!(a::InvalidState4, b::InvalidState4) = copyto!(a.Ψ, b.Ψ)
    LinearAlgebra.dot(a::InvalidState4, b::InvalidState4) = a.Ψ ⋅ b.Ψ
    LinearAlgebra.norm(a::InvalidState4) = norm(a.Ψ)
    Base.:+(a::InvalidState4, b::InvalidState4) = InvalidState4(a.Ψ + b.Ψ)
    Base.:-(a::InvalidState4, b::InvalidState4) = InvalidState4(a.Ψ - b.Ψ)
    Base.:*(α::Number, state::InvalidState4) = InvalidState4(α * state.Ψ)
    state = InvalidState4([1im, 0])
    captured = IOCapture.capture() do
        check_state(state; normalized=true)
    end
    @test captured.value ≡ false
    @test contains(captured.output, "`lmul!(c, state)` for a scalar `c` must be defined")
    @test contains(captured.output, "`lmul!(0.0, state)` must produce a state with norm 0")
    @test contains(captured.output, "`axpy!(c, state, other)` must be defined")

    state = [1, 0, 0, 0]
    captured = IOCapture.capture() do
        check_state(state; normalized=true)
    end
    @test captured.value ≡ false
    @test contains(captured.output, "`state ⋅ state` must return a Complex number type")

end


@testset "Invalid propagator" begin

    struct InvalidPropagatorEmpty end

    include("invalid_propagator.jl")

    propagator = InvalidPropagatorEmpty()

    captured = IOCapture.capture() do
        check_propagator(propagator)
    end
    @test captured.value ≡ false

    @test contains(captured.output, "does not have the required properties")

    N = 10
    tlist = collect(range(0, 100, length=1001))
    rng = StableRNG(93655235)
    Ĥ = random_dynamic_generator(N, tlist; rng)
    Ψ = random_state_vector(N; rng)

    propagator = init_prop(Ψ, Ĥ, tlist; method=:invalid_propagator_no_methods)

    captured = IOCapture.capture() do
        check_propagator(propagator)
    end

    @test contains(captured.output, "`prop_step!(propagator)` must be defined")
    @test contains(captured.output, "Failed to run `prop_step!(propagator)`")
    @test contains(captured.output, "`reinit_prop!` must be defined")

    propagator = init_prop(
        Ψ,
        Ĥ,
        tlist;
        method=:invalid_propagator_empty_methods,
        inplace=true,
        backward=false
    )
    captured = IOCapture.capture() do
        check_propagator(propagator)
    end
    @test captured.value ≡ false
    @test contains(captured.output, "propagator.t ≠ propagator.tlist[begin]")
    @test contains(captured.output, "`set_t!(propagator, t)` must set propagator.t")
    @test contains(
        captured.output,
        "`prop_step!(propagator)` at final t=0.1 must return `nothing`"
    )
    @test contains(captured.output, "`propagator.parameters` must be a dict")
    @test contains(
        captured.output,
        "`reinit_prop!(propagator, state)` must reset `propagator.t`"
    )

    propagator = init_prop(
        Ψ,
        Ĥ,
        tlist;
        method=:invalid_propagator_empty_methods,
        inplace=false,
        backward=true
    )
    captured = IOCapture.capture() do
        check_propagator(propagator)
    end
    @test captured.value ≡ false
    @test contains(captured.output, "propagator.t ≠ propagator.tlist[end]")
    @test contains(
        captured.output,
        "For a not-in-place propagator, the state returned by `prop_step!` must be a new object"
    )
    @test contains(
        captured.output,
        "`prop_step!` must advance `propagator.t` forward or backward one step on the time grid"
    )

    propagator = init_prop(
        Ψ,
        Ĥ,
        tlist;
        method=:invalid_random_propagator,
        inplace=true,
        backward=false
    )
    captured = IOCapture.capture() do
        check_propagator(propagator)
    end
    @test captured.value ≡ false
    @test contains(
        captured.output,
        "For an in-place propagator, the state returned by `prop_step!` must be the `propagator.state` object"
    )
    @test contains(
        captured.output,
        "`set_state!(propagator, state)` for an in-place propagator must overwrite `propagator.state` in-place."
    )
    @test contains(captured.output, "`reinit_prop!(propagator, state)` must be idempotent")

    propagator = init_prop(
        Ψ,
        Ĥ,
        tlist;
        method=:invalid_random_propagator,
        inplace=false,
        backward=true
    )
    captured = IOCapture.capture() do
        check_propagator(propagator)
    end
    @test captured.value ≡ false
    @test contains(
        captured.output,
        "For a not-in-place propagator, the state returned by `prop_step!` must be a new object"
    )
    @test contains(
        captured.output,
        "`reinit_prop!` must be defined and re-initialize the propagator"
    )

    propagator = init_prop(
        Ψ,
        Ĥ,
        tlist;
        method=:invalid_propagator_no_state,
        inplace=true,
        backward=false
    )
    captured = IOCapture.capture() do
        check_propagator(propagator)
    end
    @test captured.value ≡ false
    @test contains(
        captured.output,
        "prop_step! must return a valid state until time grid is exhausted"
    )

end
