using Test
using Logging: with_logger
using QuantumControlTestUtils: QuantumTestLogger
using QuantumPropagators.Interfaces:
    check_amplitude, check_control, check_operator, check_generator, check_state


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
    @test "`substitute(ampl, replacments)` must be defined" ∈ test_logger
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
    @test "`copyto!(other, state)` must be defined" ∈ test_logger
    @test "`lmul!(c, state)` for a scalar `c` must be defined" ∈ test_logger
    @test "`lmul!(0.0, state)` must produce a state with norm 0" ∈ test_logger
    @test "`axpy!(c, state, other)` must be defined" ∈ test_logger
    @test "`norm(state)` must be 1" ∈ test_logger


    state = [1, 0, 0, 0]
    test_logger = QuantumTestLogger()
    with_logger(test_logger) do
        @test check_state(state; normalized=true) ≡ false
    end
    @test "`state ⋅ state` must return a Complex number type" ∈ test_logger

end
