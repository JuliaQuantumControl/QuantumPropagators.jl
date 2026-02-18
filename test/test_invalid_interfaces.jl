using Test
using Logging: with_logger
using IOCapture: IOCapture
using QuantumControlTestUtils.RandomObjects: random_dynamic_generator, random_state_vector
using StableRNGs: StableRNG
using QuantumPropagators: QuantumPropagators, init_prop
using LinearAlgebra
import QuantumPropagators.Interfaces: supports_inplace
using QuantumPropagators.Interfaces:
    supports_matrix_interface,
    supports_vector_interface,
    check_amplitude,
    check_control,
    check_operator,
    check_generator,
    check_state,
    check_propagator

# Helper type: an immutable object with wrong size/eltype, used by "bad similar" tests
struct ImmutableResult end
QuantumPropagators.Interfaces.supports_inplace(::Type{ImmutableResult}) = false
Base.size(::ImmutableResult) = (1,)
Base.eltype(::Type{ImmutableResult}) = Float32


@testset "Invalid amplitude" begin

    struct InvalidAmplitude end

    tlist = collect(range(0, 10, length = 101))

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

    tlist = collect(range(0, 10, length = 101))

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
    tlist = collect(range(0, 10, length = 101))

    operator = InvalidOperator()
    captured = IOCapture.capture() do
        check_operator(operator; state, tlist)
    end
    @test captured.value ≡ false
    @test contains(
        captured.output,
        "The `QuantumPropagators.Interfaces.supports_inplace` method must be defined for type"
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


@testset "Invalid operator with matrix interface" begin

    struct InvalidMatrixOp end
    QuantumPropagators.Interfaces.supports_inplace(::Type{InvalidMatrixOp}) = true
    QuantumPropagators.Interfaces.supports_matrix_interface(::Type{InvalidMatrixOp}) = true
    Base.size(::InvalidMatrixOp) = (4, 4)
    Base.size(::InvalidMatrixOp, d::Int) = 4
    QuantumPropagators.Controls.get_controls(::InvalidMatrixOp) = ()
    QuantumPropagators.Controls.evaluate(op::InvalidMatrixOp, args...; kwargs...) = op
    Base.:(*)(::InvalidMatrixOp, Ψ::Vector{ComplexF64}) = Ψ
    LinearAlgebra.mul!(ϕ, ::InvalidMatrixOp, Ψ) = copyto!(ϕ, Ψ)
    LinearAlgebra.mul!(ϕ, ::InvalidMatrixOp, Ψ, α, β) = (lmul!(β, ϕ); axpy!(α, Ψ, ϕ))
    LinearAlgebra.dot(a, ::InvalidMatrixOp, b) = dot(a, b)

    state = ComplexF64[1, 0, 0, 0]
    tlist = collect(range(0, 10, length = 101))
    operator = InvalidMatrixOp()

    captured = IOCapture.capture() do
        check_operator(operator; state, tlist)
    end
    @test captured.value ≡ false
    @test contains(captured.output, "`eltype(op)` must return a numeric type")
    @test contains(captured.output, "`getindex(op, i, j)` must be defined")
    @test contains(captured.output, "`length(op)` must be defined")
    @test contains(captured.output, "`iterate(op)` must be defined")
    @test contains(captured.output, "`similar(op)` must be defined")
    @test contains(captured.output, "`similar(op, ::Type{S})` must be defined")
    @test contains(captured.output, "`similar(op, dims::Dims)` must be defined")
    @test contains(captured.output, "`similar(op, ::Type{S}, dims::Dims)` must be defined")

end


@testset "Invalid operator with bad similar (matrix interface)" begin

    # Operator that fully implements the matrix interface, but `similar` returns
    # an immutable object with wrong size and eltype.
    struct BadSimilarMatrixOp
        data::Matrix{ComplexF64}
    end
    QuantumPropagators.Interfaces.supports_inplace(::Type{BadSimilarMatrixOp}) = true
    QuantumPropagators.Interfaces.supports_matrix_interface(::Type{BadSimilarMatrixOp}) =
        true
    Base.size(op::BadSimilarMatrixOp) = size(op.data)
    Base.size(op::BadSimilarMatrixOp, d::Int) = size(op.data, d)
    Base.eltype(::Type{BadSimilarMatrixOp}) = ComplexF64
    Base.getindex(op::BadSimilarMatrixOp, i, j) = op.data[i, j]
    Base.length(op::BadSimilarMatrixOp) = length(op.data)
    Base.iterate(op::BadSimilarMatrixOp, args...) = iterate(op.data, args...)
    Base.similar(::BadSimilarMatrixOp, args...) = ImmutableResult()
    QuantumPropagators.Controls.get_controls(::BadSimilarMatrixOp) = ()
    QuantumPropagators.Controls.evaluate(op::BadSimilarMatrixOp, args...; kwargs...) = op
    Base.:(*)(op::BadSimilarMatrixOp, Ψ) = op.data * Ψ
    LinearAlgebra.mul!(ϕ, op::BadSimilarMatrixOp, Ψ) = mul!(ϕ, op.data, Ψ)
    LinearAlgebra.mul!(ϕ, op::BadSimilarMatrixOp, Ψ, α, β) = mul!(ϕ, op.data, Ψ, α, β)
    LinearAlgebra.dot(a, op::BadSimilarMatrixOp, b) = dot(a, op.data, b)

    state = ComplexF64[1, 0, 0, 0]
    tlist = collect(range(0, 10, length = 101))
    H = zeros(ComplexF64, 4, 4)
    H[1, 1] = 1.0
    operator = BadSimilarMatrixOp(H)

    captured = IOCapture.capture() do
        check_operator(operator; state, tlist)
    end
    @test captured.value ≡ false
    # similar(op): wrong mutability, shape, and eltype
    @test contains(
        captured.output,
        "`similar(op)` must return a mutable array (`supports_inplace` must be `true`)"
    )
    @test contains(
        captured.output,
        "`similar(op)` must return an array with the same shape"
    )
    @test contains(
        captured.output,
        "`similar(op)` must return an array with the same element type"
    )
    # similar(op, S):
    @test contains(captured.output, "`similar(op, ComplexF32)` must return a mutable array")
    @test contains(
        captured.output,
        "`similar(op, ComplexF32)` must return an array with the same shape"
    )
    @test contains(
        captured.output,
        "`similar(op, ComplexF32)` must return an array with element type ComplexF32"
    )
    # similar(op, dims):
    @test contains(captured.output, "`similar(op, dims)` must return a mutable array")
    @test contains(captured.output, "`similar(op, dims)` must return an array with size")
    @test contains(
        captured.output,
        "`similar(op, dims)` must return an array with the same element type"
    )
    # similar(op, S, dims):
    @test contains(
        captured.output,
        "`similar(op, ComplexF32, dims)` must return a mutable array"
    )
    @test contains(
        captured.output,
        "`similar(op, ComplexF32, dims)` must return an array with size"
    )
    @test contains(
        captured.output,
        "`similar(op, ComplexF32, dims)` must return an array with element type ComplexF32"
    )

end


@testset "Invalid operator with wrong returns" begin

    # Operator where `*` returns wrong type, `mul!` returns wrong reference,
    # and `dot` returns a mismatched value.
    struct WrongReturnOp end
    QuantumPropagators.Interfaces.supports_inplace(::Type{WrongReturnOp}) = true
    QuantumPropagators.Interfaces.supports_matrix_interface(::Type{WrongReturnOp}) = false
    Base.size(::WrongReturnOp) = (4, 4)
    Base.size(::WrongReturnOp, d::Int) = 4
    QuantumPropagators.Controls.get_controls(::WrongReturnOp) = ()
    QuantumPropagators.Controls.evaluate(op::WrongReturnOp, args...; kwargs...) = op
    Base.:(*)(::WrongReturnOp, Ψ::Vector{ComplexF64}) = real.(Ψ)  # wrong type
    function LinearAlgebra.mul!(ϕ, ::WrongReturnOp, Ψ)
        ϕ .= 3 .* Ψ
        return Ψ  # wrong: should return ϕ
    end
    function LinearAlgebra.mul!(ϕ, ::WrongReturnOp, Ψ, α, β)
        ϕ .= 3 .* Ψ
        return Ψ  # wrong: should return ϕ
    end
    LinearAlgebra.dot(a, ::WrongReturnOp, b) = 999.0 + 0im  # wrong value

    state = ComplexF64[1, 0, 0, 0]
    tlist = collect(range(0, 10, length = 101))
    operator = WrongReturnOp()

    captured = IOCapture.capture() do
        check_operator(operator; state, tlist)
    end
    @test captured.value ≡ false
    @test contains(
        captured.output,
        "`op * state` must return an object of the same type as `state`"
    )
    @test contains(captured.output, "`mul!(ϕ, op, state)` must return the resulting ϕ")
    @test contains(captured.output, "`mul!(ϕ, op, state)` must match `op * state`")
    @test contains(
        captured.output,
        "`mul!(ϕ, op, state, α, β)` must return the resulting ϕ"
    )
    @test contains(
        captured.output,
        "`mul!(ϕ, op, state, α, β)` must match β*ϕ + α*op*state"
    )
    @test contains(
        captured.output,
        "`dot(state, op, state)` must match `dot(state, op * state)`"
    )

end


@testset "Invalid operator with bad size" begin

    struct BadSizeOp end
    QuantumPropagators.Interfaces.supports_inplace(::Type{BadSizeOp}) = true
    Base.size(::BadSizeOp) = 4  # not a Tuple!
    Base.size(::BadSizeOp, d::Int) = 4
    QuantumPropagators.Controls.get_controls(::BadSizeOp) = ()
    QuantumPropagators.Controls.evaluate(op::BadSizeOp, args...; kwargs...) = op
    Base.:(*)(::BadSizeOp, Ψ::Vector{ComplexF64}) = Ψ
    LinearAlgebra.mul!(ϕ, ::BadSizeOp, Ψ) = copyto!(ϕ, Ψ)
    LinearAlgebra.mul!(ϕ, ::BadSizeOp, Ψ, α, β) = (lmul!(β, ϕ); axpy!(α, Ψ, ϕ))
    LinearAlgebra.dot(a, ::BadSizeOp, b) = dot(a, b)

    state = ComplexF64[1, 0, 0, 0]
    tlist = collect(range(0, 10, length = 101))

    captured = IOCapture.capture() do
        check_operator(BadSizeOp(); state, tlist)
    end
    @test captured.value ≡ false
    @test contains(captured.output, "`size(op)` must return a tuple")

    struct BadSizeOp2 end
    QuantumPropagators.Interfaces.supports_inplace(::Type{BadSizeOp2}) = true
    Base.size(::BadSizeOp2) = (1.5, 2.5)  # non-integer tuple
    Base.size(::BadSizeOp2, d::Int) = [1.5, 2.5][d]
    QuantumPropagators.Controls.get_controls(::BadSizeOp2) = ()
    QuantumPropagators.Controls.evaluate(op::BadSizeOp2, args...; kwargs...) = op
    Base.:(*)(::BadSizeOp2, Ψ::Vector{ComplexF64}) = Ψ
    LinearAlgebra.mul!(ϕ, ::BadSizeOp2, Ψ) = copyto!(ϕ, Ψ)
    LinearAlgebra.mul!(ϕ, ::BadSizeOp2, Ψ, α, β) = (lmul!(β, ϕ); axpy!(α, Ψ, ϕ))
    LinearAlgebra.dot(a, ::BadSizeOp2, b) = dot(a, b)

    captured = IOCapture.capture() do
        check_operator(BadSizeOp2(); state, tlist)
    end
    @test captured.value ≡ false
    @test contains(captured.output, "`size(op)` must return a tuple of integers")

end


@testset "Invalid operator with bad size dimensions" begin

    # size(op, dim) returns non-integer
    struct BadSizeDimOp end
    QuantumPropagators.Interfaces.supports_inplace(::Type{BadSizeDimOp}) = true
    Base.size(::BadSizeDimOp) = (4, 4)
    Base.size(::BadSizeDimOp, d::Int) = 4.0  # Float64, not Int
    QuantumPropagators.Controls.get_controls(::BadSizeDimOp) = ()
    QuantumPropagators.Controls.evaluate(op::BadSizeDimOp, args...; kwargs...) = op
    Base.:(*)(::BadSizeDimOp, Ψ::Vector{ComplexF64}) = Ψ
    LinearAlgebra.mul!(ϕ, ::BadSizeDimOp, Ψ) = copyto!(ϕ, Ψ)
    LinearAlgebra.mul!(ϕ, ::BadSizeDimOp, Ψ, α, β) = (lmul!(β, ϕ); axpy!(α, Ψ, ϕ))
    LinearAlgebra.dot(a, ::BadSizeDimOp, b) = dot(a, b)

    state = ComplexF64[1, 0, 0, 0]
    tlist = collect(range(0, 10, length = 101))

    captured = IOCapture.capture() do
        check_operator(BadSizeDimOp(); state, tlist)
    end
    @test captured.value ≡ false
    @test contains(captured.output, "`size(op, 1)` must return an integer, not Float64")

    # size(op, dim) inconsistent with size(op)
    struct BadSizeDimOp2 end
    QuantumPropagators.Interfaces.supports_inplace(::Type{BadSizeDimOp2}) = true
    Base.size(::BadSizeDimOp2) = (4, 4)
    Base.size(::BadSizeDimOp2, d::Int) = 5  # inconsistent: 5 ≠ 4
    QuantumPropagators.Controls.get_controls(::BadSizeDimOp2) = ()
    QuantumPropagators.Controls.evaluate(op::BadSizeDimOp2, args...; kwargs...) = op
    Base.:(*)(::BadSizeDimOp2, Ψ::Vector{ComplexF64}) = Ψ
    LinearAlgebra.mul!(ϕ, ::BadSizeDimOp2, Ψ) = copyto!(ϕ, Ψ)
    LinearAlgebra.mul!(ϕ, ::BadSizeDimOp2, Ψ, α, β) = (lmul!(β, ϕ); axpy!(α, Ψ, ϕ))
    LinearAlgebra.dot(a, ::BadSizeDimOp2, b) = dot(a, b)

    captured = IOCapture.capture() do
        check_operator(BadSizeDimOp2(); state, tlist)
    end
    @test captured.value ≡ false
    @test contains(captured.output, "`size(op, 1)` must be consistent with `size(op)`")

end


@testset "Invalid operator with partial matrix interface" begin

    # Operator that has matrix interface but getindex returns wrong type,
    # length is wrong, iterate returns nothing, and get_controls is non-empty.
    struct PartialMatrixOp end
    QuantumPropagators.Interfaces.supports_inplace(::Type{PartialMatrixOp}) = false
    QuantumPropagators.Interfaces.supports_matrix_interface(::Type{PartialMatrixOp}) = true
    Base.size(::PartialMatrixOp) = (4, 4)
    Base.size(::PartialMatrixOp, d::Int) = 4
    Base.eltype(::Type{PartialMatrixOp}) = ComplexF64
    Base.getindex(::PartialMatrixOp, i, j) = 1.0  # Float64, not ComplexF64!
    Base.length(::PartialMatrixOp) = 99  # wrong: should be 16
    Base.iterate(::PartialMatrixOp) = nothing  # wrong for non-empty
    Base.similar(::PartialMatrixOp, args...) = zeros(ComplexF64, 4, 4)  # correct
    QuantumPropagators.Controls.get_controls(::PartialMatrixOp) = (1.0,)  # non-empty!
    QuantumPropagators.Controls.evaluate(op::PartialMatrixOp, args...; kwargs...) = op
    Base.:(*)(::PartialMatrixOp, Ψ::Vector{ComplexF64}) = Ψ
    LinearAlgebra.dot(a, ::PartialMatrixOp, b) = dot(a, b)

    state = ComplexF64[1, 0, 0, 0]
    tlist = collect(range(0, 10, length = 101))

    captured = IOCapture.capture() do
        check_operator(PartialMatrixOp(); state, tlist)
    end
    @test captured.value ≡ false
    @test contains(captured.output, "get_controls(op) must return an empty tuple")
    @test contains(
        captured.output,
        "`op[1, 1]` must return a value of type `eltype(op)=ComplexF64`"
    )
    @test contains(captured.output, "`length(op)` must equal `prod(size(op))`")
    @test contains(
        captured.output,
        "`iterate(op)` must not return `nothing` for a non-empty operator"
    )

end


@testset "Invalid generator" begin

    struct InvalidGenerator end

    state = ComplexF64[1, 0, 0, 0]
    tlist = collect(range(0, 10, length = 101))

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
        "The `QuantumPropagators.Interfaces.supports_inplace` method must be defined for type"
    )

    struct InvalidState end
    QuantumPropagators.Interfaces.supports_inplace(::Type{InvalidState}) = true
    state = InvalidState()
    captured = IOCapture.capture() do
        check_state(state; normalized = true)
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
    QuantumPropagators.Interfaces.supports_inplace(::Type{InvalidState2}) = true
    state = InvalidState2()
    Base.similar(::InvalidState2) = InvalidState2()
    captured = IOCapture.capture() do
        check_state(state; normalized = true)
    end
    @test captured.value ≡ false
    @test contains(captured.output, "`copyto!(other, state)` must be defined")

    struct InvalidState3 end
    QuantumPropagators.Interfaces.supports_inplace(::Type{InvalidState3}) = true
    state = InvalidState3()
    Base.similar(::InvalidState3) = InvalidState3()
    Base.copyto!(a::InvalidState3, b::InvalidState3) = a
    captured = IOCapture.capture() do
        check_state(state; normalized = true)
    end
    @test captured.value ≡ false
    @test contains(captured.output, "`similar(state)` must return a valid state")
    @test contains(captured.output, "On `similar(state)`: ")

    struct InvalidState4
        Ψ
    end
    QuantumPropagators.Interfaces.supports_inplace(::Type{InvalidState4}) = true
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
        check_state(state; normalized = true)
    end
    @test captured.value ≡ false
    @test contains(captured.output, "`lmul!(c, state)` for a scalar `c` must be defined")
    @test contains(captured.output, "`lmul!(0.0, state)` must produce a state with norm 0")
    @test contains(captured.output, "`axpy!(c, state, other)` must be defined")

    state = [1, 0, 0, 0]
    captured = IOCapture.capture() do
        check_state(state; normalized = true)
    end
    @test captured.value ≡ false
    @test contains(captured.output, "`state ⋅ state` must return a Complex number type")

end


@testset "Invalid state with vector interface" begin

    struct InvalidVectorState
        data::Vector{ComplexF64}
    end
    QuantumPropagators.Interfaces.supports_inplace(::Type{InvalidVectorState}) = true
    QuantumPropagators.Interfaces.supports_vector_interface(::Type{InvalidVectorState}) =
        true
    Base.copy(s::InvalidVectorState) = InvalidVectorState(copy(s.data))
    Base.similar(s::InvalidVectorState) = InvalidVectorState(similar(s.data))
    Base.copyto!(a::InvalidVectorState, b::InvalidVectorState) =
        (copyto!(a.data, b.data); a)
    LinearAlgebra.dot(a::InvalidVectorState, b::InvalidVectorState) = dot(a.data, b.data)
    LinearAlgebra.norm(a::InvalidVectorState) = norm(a.data)
    Base.:+(a::InvalidVectorState, b::InvalidVectorState) =
        InvalidVectorState(a.data + b.data)
    Base.:-(a::InvalidVectorState, b::InvalidVectorState) =
        InvalidVectorState(a.data - b.data)
    Base.:*(α::Number, s::InvalidVectorState) = InvalidVectorState(α * s.data)
    Base.zero(s::InvalidVectorState) = InvalidVectorState(zero(s.data))
    Base.fill!(s::InvalidVectorState, v) = (fill!(s.data, v); s)
    LinearAlgebra.lmul!(c, s::InvalidVectorState) = (lmul!(c, s.data); s)
    LinearAlgebra.axpy!(c, a::InvalidVectorState, b::InvalidVectorState) =
        (axpy!(c, a.data, b.data); b)
    # Deliberately don't define: eltype, getindex, setindex!, length, iterate, size

    state = InvalidVectorState(ComplexF64[1, 0, 0, 0])
    captured = IOCapture.capture() do
        check_state(state; normalized = true)
    end
    @test captured.value ≡ false
    @test contains(captured.output, "`eltype(state)` must return a numeric type")
    @test contains(captured.output, "`getindex(state, i)` must be defined")
    @test contains(captured.output, "`length(state)` must be defined")
    @test contains(captured.output, "`iterate(state)` must be defined")
    @test contains(captured.output, "`setindex!(state, v, i)` must be defined")

end


@testset "Invalid state with vector interface (no similar)" begin

    struct InvalidVectorState2
        data::Vector{ComplexF64}
    end
    QuantumPropagators.Interfaces.supports_inplace(::Type{InvalidVectorState2}) = false
    QuantumPropagators.Interfaces.supports_vector_interface(::Type{InvalidVectorState2}) =
        true
    Base.copy(s::InvalidVectorState2) = InvalidVectorState2(copy(s.data))
    LinearAlgebra.dot(a::InvalidVectorState2, b::InvalidVectorState2) = dot(a.data, b.data)
    LinearAlgebra.norm(a::InvalidVectorState2) = norm(a.data)
    Base.:+(a::InvalidVectorState2, b::InvalidVectorState2) =
        InvalidVectorState2(a.data + b.data)
    Base.:-(a::InvalidVectorState2, b::InvalidVectorState2) =
        InvalidVectorState2(a.data - b.data)
    Base.:*(α::Number, s::InvalidVectorState2) = InvalidVectorState2(α * s.data)
    Base.zero(s::InvalidVectorState2) = InvalidVectorState2(zero(s.data))
    # Deliberately don't define: eltype, getindex, length, iterate, similar, size

    state = InvalidVectorState2(ComplexF64[1, 0, 0, 0])
    captured = IOCapture.capture() do
        check_state(state; normalized = true)
    end
    @test captured.value ≡ false
    @test contains(captured.output, "`eltype(state)` must return a numeric type")
    @test contains(captured.output, "`getindex(state, i)` must be defined")
    @test contains(captured.output, "`length(state)` must be defined")
    @test contains(captured.output, "`iterate(state)` must be defined")
    @test contains(captured.output, "`similar(state)` must be defined")
    @test contains(captured.output, "`similar(state, ::Type{S})` must be defined")
    @test contains(captured.output, "`similar(state, dims::Dims)` must be defined")
    @test contains(
        captured.output,
        "`similar(state, ::Type{S}, dims::Dims)` must be defined"
    )

end


@testset "Invalid state with bad similar (vector interface)" begin

    # State that fully implements the Hilbert space and array interface, but
    # `similar` returns an immutable object with wrong size and eltype.
    struct BadSimilarVectorState3
        data::Vector{ComplexF64}
    end
    QuantumPropagators.Interfaces.supports_inplace(::Type{BadSimilarVectorState3}) = false
    QuantumPropagators.Interfaces.supports_vector_interface(
        ::Type{BadSimilarVectorState3}
    ) = true
    Base.copy(s::BadSimilarVectorState3) = BadSimilarVectorState3(copy(s.data))
    LinearAlgebra.dot(a::BadSimilarVectorState3, b::BadSimilarVectorState3) =
        dot(a.data, b.data)
    LinearAlgebra.norm(a::BadSimilarVectorState3) = norm(a.data)
    Base.:+(a::BadSimilarVectorState3, b::BadSimilarVectorState3) =
        BadSimilarVectorState3(a.data + b.data)
    Base.:-(a::BadSimilarVectorState3, b::BadSimilarVectorState3) =
        BadSimilarVectorState3(a.data - b.data)
    Base.:*(α::Number, s::BadSimilarVectorState3) = BadSimilarVectorState3(α * s.data)
    Base.zero(s::BadSimilarVectorState3) = BadSimilarVectorState3(zero(s.data))
    Base.eltype(::Type{BadSimilarVectorState3}) = ComplexF64
    Base.size(s::BadSimilarVectorState3) = size(s.data)
    Base.getindex(s::BadSimilarVectorState3, i) = s.data[i]
    Base.length(s::BadSimilarVectorState3) = length(s.data)
    Base.iterate(s::BadSimilarVectorState3, args...) = iterate(s.data, args...)
    Base.similar(::BadSimilarVectorState3, args...) = ImmutableResult()

    state = BadSimilarVectorState3(ComplexF64[1, 0, 0, 0])
    captured = IOCapture.capture() do
        check_state(state; normalized = true)
    end
    @test captured.value ≡ false
    # similar(state): wrong mutability, shape, and eltype
    @test contains(
        captured.output,
        "`similar(state)` must return a mutable vector (`supports_inplace` must be `true`)"
    )
    @test contains(
        captured.output,
        "`similar(state)` must return a vector with the same shape"
    )
    @test contains(
        captured.output,
        "`similar(state)` must return a vector with the same element type"
    )
    # similar(state, S):
    @test contains(
        captured.output,
        "`similar(state, ComplexF32)` must return a mutable vector"
    )
    @test contains(
        captured.output,
        "`similar(state, ComplexF32)` must return a vector with the same shape"
    )
    @test contains(
        captured.output,
        "`similar(state, ComplexF32)` must return a vector with element type ComplexF32"
    )
    # similar(state, dims):
    @test contains(captured.output, "`similar(state, dims)` must return a mutable array")
    @test contains(captured.output, "`similar(state, dims)` must return an array with size")
    @test contains(
        captured.output,
        "`similar(state, dims)` must return an array with the same element type"
    )
    # similar(state, S, dims):
    @test contains(
        captured.output,
        "`similar(state, ComplexF32, dims)` must return a mutable array"
    )
    @test contains(
        captured.output,
        "`similar(state, ComplexF32, dims)` must return an array with size"
    )
    @test contains(
        captured.output,
        "`similar(state, ComplexF32, dims)` must return an array with element type"
    )

end


@testset "Invalid state non-unit norm" begin

    state = ComplexF64[2, 0, 0, 0]
    captured = IOCapture.capture() do
        check_state(state; normalized = true)
    end
    @test captured.value ≡ false
    @test contains(captured.output, "`norm(state)` must be 1, not 2.0")

end


@testset "Invalid state with partial vector interface" begin

    # State that supports vector interface but getindex returns wrong type,
    # length is wrong, and iterate returns nothing.
    struct PartialVectorState4
        data::Vector{ComplexF64}
    end
    QuantumPropagators.Interfaces.supports_inplace(::Type{PartialVectorState4}) = false
    QuantumPropagators.Interfaces.supports_vector_interface(::Type{PartialVectorState4}) =
        true
    Base.copy(s::PartialVectorState4) = PartialVectorState4(copy(s.data))
    LinearAlgebra.dot(a::PartialVectorState4, b::PartialVectorState4) = dot(a.data, b.data)
    LinearAlgebra.norm(a::PartialVectorState4) = norm(a.data)
    Base.:+(a::PartialVectorState4, b::PartialVectorState4) =
        PartialVectorState4(a.data + b.data)
    Base.:-(a::PartialVectorState4, b::PartialVectorState4) =
        PartialVectorState4(a.data - b.data)
    Base.:*(α::Number, s::PartialVectorState4) = PartialVectorState4(α * s.data)
    Base.zero(s::PartialVectorState4) = PartialVectorState4(zero(s.data))
    Base.eltype(::Type{PartialVectorState4}) = ComplexF64
    Base.size(s::PartialVectorState4) = size(s.data)
    Base.getindex(s::PartialVectorState4, i) = real(s.data[i])  # Float64, not ComplexF64!
    Base.length(::PartialVectorState4) = 99  # wrong: should be 4
    Base.iterate(::PartialVectorState4) = nothing  # wrong for non-empty
    Base.similar(::PartialVectorState4, args...) = zeros(ComplexF64, 4)

    state = PartialVectorState4(ComplexF64[1, 0, 0, 0])
    captured = IOCapture.capture() do
        check_state(state; normalized = true)
    end
    @test captured.value ≡ false
    @test contains(
        captured.output,
        "`state[1]` must return a value of type `eltype(state)=ComplexF64`"
    )
    @test contains(captured.output, "`length(state)` must equal `prod(size(state))`")
    @test contains(
        captured.output,
        "`iterate(state)` must not return `nothing` for a non-empty state"
    )

end


@testset "Invalid operator with throwing size" begin

    # Operator where size(op) itself throws
    struct ThrowingSizeOp end
    QuantumPropagators.Interfaces.supports_inplace(::Type{ThrowingSizeOp}) = true
    Base.size(::ThrowingSizeOp) = error("size not implemented")
    QuantumPropagators.Controls.get_controls(::ThrowingSizeOp) = ()
    QuantumPropagators.Controls.evaluate(op::ThrowingSizeOp, args...; kwargs...) = op
    Base.:(*)(::ThrowingSizeOp, Ψ::Vector{ComplexF64}) = Ψ
    LinearAlgebra.mul!(ϕ, ::ThrowingSizeOp, Ψ) = copyto!(ϕ, Ψ)
    LinearAlgebra.mul!(ϕ, ::ThrowingSizeOp, Ψ, α, β) = (lmul!(β, ϕ); axpy!(α, Ψ, ϕ))
    LinearAlgebra.dot(a, ::ThrowingSizeOp, b) = dot(a, b)

    state = ComplexF64[1, 0, 0, 0]
    tlist = collect(range(0, 10, length = 101))

    captured = IOCapture.capture() do
        check_operator(ThrowingSizeOp(); state, tlist)
    end
    @test captured.value ≡ false
    @test contains(captured.output, "`size(op)` must be defined")

end


@testset "Invalid operator with throwing evaluate" begin

    # Operator where evaluate throws
    struct ThrowingEvalOp end
    QuantumPropagators.Interfaces.supports_inplace(::Type{ThrowingEvalOp}) = true
    Base.size(::ThrowingEvalOp) = (4, 4)
    Base.size(::ThrowingEvalOp, d::Int) = 4
    QuantumPropagators.Controls.get_controls(::ThrowingEvalOp) = ()
    QuantumPropagators.Controls.evaluate(::ThrowingEvalOp, args...; kwargs...) =
        error("cannot evaluate")
    Base.:(*)(::ThrowingEvalOp, Ψ::Vector{ComplexF64}) = Ψ
    LinearAlgebra.mul!(ϕ, ::ThrowingEvalOp, Ψ) = copyto!(ϕ, Ψ)
    LinearAlgebra.mul!(ϕ, ::ThrowingEvalOp, Ψ, α, β) = (lmul!(β, ϕ); axpy!(α, Ψ, ϕ))
    LinearAlgebra.dot(a, ::ThrowingEvalOp, b) = dot(a, b)

    state = ComplexF64[1, 0, 0, 0]
    tlist = collect(range(0, 10, length = 101))

    captured = IOCapture.capture() do
        check_operator(ThrowingEvalOp(); state, tlist)
    end
    @test captured.value ≡ false
    @test contains(captured.output, "op must not be time-dependent")

end


@testset "Invalid operator with throwing size(op, dim)" begin

    # Operator where size(op) works but size(op, dim) throws
    struct ThrowingSizeDimOp end
    QuantumPropagators.Interfaces.supports_inplace(::Type{ThrowingSizeDimOp}) = true
    Base.size(::ThrowingSizeDimOp) = (4, 4)
    Base.size(::ThrowingSizeDimOp, d::Int) = error("size not implemented for dim")
    QuantumPropagators.Controls.get_controls(::ThrowingSizeDimOp) = ()
    QuantumPropagators.Controls.evaluate(op::ThrowingSizeDimOp, args...; kwargs...) = op
    Base.:(*)(::ThrowingSizeDimOp, Ψ::Vector{ComplexF64}) = Ψ
    LinearAlgebra.mul!(ϕ, ::ThrowingSizeDimOp, Ψ) = copyto!(ϕ, Ψ)
    LinearAlgebra.mul!(ϕ, ::ThrowingSizeDimOp, Ψ, α, β) = (lmul!(β, ϕ); axpy!(α, Ψ, ϕ))
    LinearAlgebra.dot(a, ::ThrowingSizeDimOp, b) = dot(a, b)

    state = ComplexF64[1, 0, 0, 0]
    tlist = collect(range(0, 10, length = 101))

    captured = IOCapture.capture() do
        check_operator(ThrowingSizeDimOp(); state, tlist)
    end
    @test captured.value ≡ false
    @test contains(captured.output, "`size(op, 1)` must be defined")

end


@testset "Invalid state with wrong copy type" begin

    # State where copy returns the wrong type (inner Vector instead of wrapper)
    struct WrongCopyState
        data::Vector{ComplexF64}
    end
    QuantumPropagators.Interfaces.supports_inplace(::Type{WrongCopyState}) = true
    Base.copy(s::WrongCopyState) = copy(s.data)  # wrong: returns Vector, not WrongCopyState
    Base.similar(s::WrongCopyState) = WrongCopyState(similar(s.data))
    Base.copyto!(a::WrongCopyState, b::WrongCopyState) = (copyto!(a.data, b.data); a)
    LinearAlgebra.dot(a::WrongCopyState, b::WrongCopyState) = dot(a.data, b.data)
    LinearAlgebra.norm(a::WrongCopyState) = norm(a.data)
    Base.:+(a::WrongCopyState, b::WrongCopyState) = WrongCopyState(a.data + b.data)
    Base.:-(a::WrongCopyState, b::WrongCopyState) = WrongCopyState(a.data - b.data)
    Base.:*(α::Number, s::WrongCopyState) = WrongCopyState(α * s.data)
    Base.zero(s::WrongCopyState) = WrongCopyState(zero(s.data))
    Base.fill!(s::WrongCopyState, v) = (fill!(s.data, v); s)
    LinearAlgebra.lmul!(c, s::WrongCopyState) = (lmul!(c, s.data); s)
    LinearAlgebra.axpy!(c, a::WrongCopyState, b::WrongCopyState) =
        (axpy!(c, a.data, b.data); b)

    state = WrongCopyState(ComplexF64[1, 0, 0, 0])
    captured = IOCapture.capture() do
        check_state(state; normalized = true)
    end
    @test captured.value ≡ false
    @test contains(
        captured.output,
        "`copy(state)::Vector{ComplexF64}` must have the same type as `state::"
    )

end


@testset "Invalid state with constant norm" begin

    # State where norm always returns a constant (inconsistent with dot)
    struct ConstantNormState
        data::Vector{ComplexF64}
    end
    QuantumPropagators.Interfaces.supports_inplace(::Type{ConstantNormState}) = false
    Base.copy(s::ConstantNormState) = ConstantNormState(copy(s.data))
    LinearAlgebra.dot(a::ConstantNormState, b::ConstantNormState) = dot(a.data, b.data)
    LinearAlgebra.norm(::ConstantNormState) = 999.0  # always 999
    Base.:+(a::ConstantNormState, b::ConstantNormState) = ConstantNormState(a.data + b.data)
    Base.:-(a::ConstantNormState, b::ConstantNormState) = ConstantNormState(a.data - b.data)
    Base.:*(α::Number, s::ConstantNormState) = ConstantNormState(α * s.data)
    Base.zero(s::ConstantNormState) = ConstantNormState(zero(s.data))

    state = ConstantNormState(ComplexF64[1, 0, 0, 0])
    captured = IOCapture.capture() do
        check_state(state)
    end
    @test captured.value ≡ false
    @test contains(captured.output, "must match `√(state⋅state)")
    @test contains(captured.output, "`state - state` must have norm 0")

end


@testset "Invalid state with squared norm" begin

    # State where norm = sum(abs2) which violates triangle inequality
    struct SquaredNormState
        data::Vector{ComplexF64}
    end
    QuantumPropagators.Interfaces.supports_inplace(::Type{SquaredNormState}) = false
    Base.copy(s::SquaredNormState) = SquaredNormState(copy(s.data))
    LinearAlgebra.dot(a::SquaredNormState, b::SquaredNormState) = dot(a.data, b.data)
    LinearAlgebra.norm(s::SquaredNormState) = sum(abs2, s.data)
    Base.:+(a::SquaredNormState, b::SquaredNormState) = SquaredNormState(a.data + b.data)
    Base.:-(a::SquaredNormState, b::SquaredNormState) = SquaredNormState(a.data - b.data)
    Base.:*(α::Number, s::SquaredNormState) = SquaredNormState(α * s.data)
    Base.zero(s::SquaredNormState) = SquaredNormState(zero(s.data))

    state = SquaredNormState(ComplexF64[1, 0, 0, 0])
    captured = IOCapture.capture() do
        check_state(state)
    end
    @test captured.value ≡ false
    @test contains(captured.output, "must fulfill the triangle inequality")

end


@testset "Invalid state with similar returning wrong type" begin

    # State where similar returns a plain Vector (wrong type) in inplace section
    struct BadSimilarTypeState
        data::Vector{ComplexF64}
    end
    QuantumPropagators.Interfaces.supports_inplace(::Type{BadSimilarTypeState}) = true
    Base.copy(s::BadSimilarTypeState) = BadSimilarTypeState(copy(s.data))
    Base.similar(s::BadSimilarTypeState) = similar(s.data)  # returns Vector, not wrapper
    Base.copyto!(a::BadSimilarTypeState, b::BadSimilarTypeState) =
        (copyto!(a.data, b.data); a)
    LinearAlgebra.dot(a::BadSimilarTypeState, b::BadSimilarTypeState) = dot(a.data, b.data)
    LinearAlgebra.norm(a::BadSimilarTypeState) = norm(a.data)
    Base.:+(a::BadSimilarTypeState, b::BadSimilarTypeState) =
        BadSimilarTypeState(a.data + b.data)
    Base.:-(a::BadSimilarTypeState, b::BadSimilarTypeState) =
        BadSimilarTypeState(a.data - b.data)
    Base.:*(α::Number, s::BadSimilarTypeState) = BadSimilarTypeState(α * s.data)
    Base.zero(s::BadSimilarTypeState) = BadSimilarTypeState(zero(s.data))
    Base.fill!(s::BadSimilarTypeState, v) = (fill!(s.data, v); s)
    LinearAlgebra.lmul!(c, s::BadSimilarTypeState) = (lmul!(c, s.data); s)
    LinearAlgebra.axpy!(c, a::BadSimilarTypeState, b::BadSimilarTypeState) =
        (axpy!(c, a.data, b.data); b)

    state = BadSimilarTypeState(ComplexF64[1, 0, 0, 0])
    captured = IOCapture.capture() do
        check_state(state; normalized = true)
    end
    @test captured.value ≡ false
    @test contains(
        captured.output,
        "`similar(state)::Vector{ComplexF64}` must have the same type as `state::"
    )

end


@testset "Invalid state with unfaithful copyto!" begin

    # State where copyto! doesn't actually copy the data
    struct BadCopytoState
        data::Vector{ComplexF64}
    end
    QuantumPropagators.Interfaces.supports_inplace(::Type{BadCopytoState}) = true
    Base.copy(s::BadCopytoState) = BadCopytoState(copy(s.data))
    Base.similar(s::BadCopytoState) = BadCopytoState(similar(s.data))
    Base.copyto!(a::BadCopytoState, b::BadCopytoState) = a  # doesn't actually copy!
    LinearAlgebra.dot(a::BadCopytoState, b::BadCopytoState) = dot(a.data, b.data)
    LinearAlgebra.norm(a::BadCopytoState) = norm(a.data)
    Base.:+(a::BadCopytoState, b::BadCopytoState) = BadCopytoState(a.data + b.data)
    Base.:-(a::BadCopytoState, b::BadCopytoState) = BadCopytoState(a.data - b.data)
    Base.:*(α::Number, s::BadCopytoState) = BadCopytoState(α * s.data)
    Base.zero(s::BadCopytoState) = BadCopytoState(zero(s.data))
    Base.fill!(s::BadCopytoState, v) = (fill!(s.data, v); s)
    LinearAlgebra.lmul!(c, s::BadCopytoState) = (lmul!(c, s.data); s)
    LinearAlgebra.axpy!(c, a::BadCopytoState, b::BadCopytoState) =
        (axpy!(c, a.data, b.data); b)

    state = BadCopytoState(ComplexF64[1, 0, 0, 0])
    captured = IOCapture.capture() do
        check_state(state; normalized = true)
    end
    @test captured.value ≡ false
    @test contains(
        captured.output,
        "must have norm 0, where `ϕ = similar(state); copyto!(ϕ, state)`"
    )

end


@testset "Invalid state with fill! not zeroing" begin

    # State where fill! returns the state (correct type) but doesn't zero the data
    struct BadFillNormState
        data::Vector{ComplexF64}
    end
    QuantumPropagators.Interfaces.supports_inplace(::Type{BadFillNormState}) = true
    Base.copy(s::BadFillNormState) = BadFillNormState(copy(s.data))
    Base.similar(s::BadFillNormState) = BadFillNormState(similar(s.data))
    Base.copyto!(a::BadFillNormState, b::BadFillNormState) = (copyto!(a.data, b.data); a)
    LinearAlgebra.dot(a::BadFillNormState, b::BadFillNormState) = dot(a.data, b.data)
    LinearAlgebra.norm(a::BadFillNormState) = norm(a.data)
    Base.:+(a::BadFillNormState, b::BadFillNormState) = BadFillNormState(a.data + b.data)
    Base.:-(a::BadFillNormState, b::BadFillNormState) = BadFillNormState(a.data - b.data)
    Base.:*(α::Number, s::BadFillNormState) = BadFillNormState(α * s.data)
    Base.zero(s::BadFillNormState) = BadFillNormState(zero(s.data))
    Base.fill!(s::BadFillNormState, v) = s  # returns state but doesn't fill!
    LinearAlgebra.lmul!(c, s::BadFillNormState) = (lmul!(c, s.data); s)
    LinearAlgebra.axpy!(c, a::BadFillNormState, b::BadFillNormState) =
        (axpy!(c, a.data, b.data); b)

    state = BadFillNormState(ComplexF64[1, 0, 0, 0])
    captured = IOCapture.capture() do
        check_state(state; normalized = true)
    end
    @test captured.value ≡ false
    @test contains(captured.output, "`fill!(state, 0.0)` must have norm 0")

end


@testset "Invalid state with buggy inplace operations" begin

    # State where fill!/lmul!/axpy! are implemented but produce wrong results
    struct BuggyInplaceState
        data::Vector{ComplexF64}
    end
    QuantumPropagators.Interfaces.supports_inplace(::Type{BuggyInplaceState}) = true
    Base.copy(s::BuggyInplaceState) = BuggyInplaceState(copy(s.data))
    Base.similar(s::BuggyInplaceState) = BuggyInplaceState(similar(s.data))
    Base.copyto!(a::BuggyInplaceState, b::BuggyInplaceState) = (copyto!(a.data, b.data); a)
    LinearAlgebra.dot(a::BuggyInplaceState, b::BuggyInplaceState) = dot(a.data, b.data)
    LinearAlgebra.norm(a::BuggyInplaceState) = norm(a.data)
    Base.:+(a::BuggyInplaceState, b::BuggyInplaceState) = BuggyInplaceState(a.data + b.data)
    Base.:-(a::BuggyInplaceState, b::BuggyInplaceState) = BuggyInplaceState(a.data - b.data)
    Base.:*(α::Number, s::BuggyInplaceState) = BuggyInplaceState(α * s.data)
    Base.zero(s::BuggyInplaceState) = BuggyInplaceState(zero(s.data))
    # fill! returns wrong type (nothing instead of the state)
    Base.fill!(s::BuggyInplaceState, v) = (fill!(s.data, v); nothing)
    # lmul! doesn't actually scale (just returns the state unchanged)
    LinearAlgebra.lmul!(c, s::BuggyInplaceState) = s
    # axpy! doesn't actually add (just returns the destination unchanged)
    LinearAlgebra.axpy!(c, a::BuggyInplaceState, b::BuggyInplaceState) = b

    state = BuggyInplaceState(ComplexF64[1, 0, 0, 0])
    captured = IOCapture.capture() do
        check_state(state; normalized = true)
    end
    @test captured.value ≡ false
    @test contains(captured.output, "`fill!(state, 0.0)` must return the filled state")
    @test contains(captured.output, "`norm(state)` must have absolute homogeneity")
    @test contains(captured.output, "`axpy!(a, state, ϕ)` must match `ϕ += a * state`")

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
    tlist = collect(range(0, 100, length = 1001))
    rng = StableRNG(93655235)
    Ĥ = random_dynamic_generator(N, tlist; rng)
    Ψ = random_state_vector(N; rng)

    propagator = init_prop(Ψ, Ĥ, tlist; method = :invalid_propagator_no_methods)

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
        method = :invalid_propagator_empty_methods,
        inplace = true,
        backward = false
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
        method = :invalid_propagator_empty_methods,
        inplace = false,
        backward = true
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
        method = :invalid_random_propagator,
        inplace = true,
        backward = false
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
        method = :invalid_random_propagator,
        inplace = false,
        backward = true
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
        method = :invalid_propagator_no_state,
        inplace = true,
        backward = false
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
