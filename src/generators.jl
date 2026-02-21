module Generators

using LinearAlgebra
using SparseArrays

export Generator, Operator, ScaledOperator
export hamiltonian, liouvillian
import ..Controls: get_controls, evaluate, evaluate!, substitute


@doc raw"""A time-dependent generator.

```julia
Generator(ops::Vector{OT}, amplitudes::Vector{AT})
```

produces an object of type `Generator{OT,AT}` that represents

```math
Ĥ(t)= Ĥ_0 + \sum_l a_l(\{ϵ_{l'}(t)\}, t) \, Ĥ_l\,,
```

where ``Ĥ_l`` are the `ops` and ``a_l(t)`` are the `amplitudes`. ``Ĥ(t)`` and
``Ĥ_l`` may represent operators in Hilbert space or super-operators in
Liouville space. If the number of `amplitudes` is less than the number of
`ops`, the first `ops` are considered as drift terms (``Ĥ_0``,
respectively subsequent terms with ``a_l ≡ 1``).
At least one time-dependent amplitude is required. Each amplitude may depend on
one or more control functions ``ϵ_{l'}(t)``, although most typically ``a_l(t) ≡
ϵ_l(t)``, that is, the `amplitudes` are simply a vector of the controls. See
[`hamiltonian`](@ref) for details.

A `Generator` object should generally not be instantiated directly, but via
[`hamiltonian`](@ref) or [`liouvillian`](@ref).

The list of `ops` and `amplitudes` are properties of the `Generator`. They
should not be mutated.

# See also

* [`Operator`](@ref) for static generators, which may be obtained from a
  `Generator` via [`evaluate`](@ref).
"""
struct Generator{OT,AT}

    ops::Vector{OT}
    amplitudes::Vector{AT}

    function Generator(ops::Vector{OT}, amplitudes::Vector{AT}) where {OT,AT}
        if length(amplitudes) > length(ops)
            error(
                "The number of amplitudes cannot exceed the number of operators in a Generator"
            )
        end
        if length(amplitudes) < 1
            error("A Generator requires at least one amplitude")
        end
        new{OT,AT}(ops, amplitudes)
    end

end

function Base.show(io::IO, G::Generator{OT,AT}) where {OT,AT}
    print(io, "Generator($(G.ops), $(G.amplitudes)")
end

function Base.summary(io::IO, G::Generator)
    print(io, "Generator with $(length(G.ops)) ops and $(length(G.amplitudes)) amplitudes")
end

function Base.show(io::IO, ::MIME"text/plain", G::Generator{OT,AT}) where {OT,AT}
    Base.summary(io, G)
    println(io, "\n ops::Vector{$OT}:")
    for op in G.ops
        print(io, "  ")
        show(io, op)
        print(io, "\n")
    end
    println(io, " amplitudes::Vector{$AT}:")
    for ampl in G.amplitudes
        print(io, "  ")
        show(io, ampl)
        print(io, "\n")
    end
end


@doc raw"""A static operator in Hilbert or Liouville space.

```julia
Operator(ops::Vector{OT}, coeffs::Vector{CT})
```

produces an object of type `Operator{OT,CT}` that encapsulates the "lazy" sum

```math
Ĥ = \sum_l c_l Ĥ_l\,,
```

where ``Ĥ_l`` are the `ops` and ``c_l`` are the `coeffs`, which each must be a
constant `Number`. If the number of coefficients is less than the
number of operators, the first `ops` are considered to have ``c_l = 1``.

An `Operator` object would generally not be instantiated directly, but be
obtained from a [`Generator`](@ref) via [`evaluate`](@ref).

The ``Ĥ_l`` in the sum are considered immutable. This implies that an
`Operator` can be updated in-place with [`evaluate!`](@ref) by only changing
the `coeffs`.
"""
struct Operator{OT,CT<:Number}

    ops::Vector{OT}
    coeffs::Vector{CT}

    function Operator(ops::Vector{OT}, coeffs::Vector{CT}) where {OT,CT}
        if length(coeffs) > length(ops)
            error(
                "The number of coefficients cannot exceed the number of operators in an Operator"
            )
        end
        new{OT,CT}(ops, coeffs)
    end

end

function Base.show(io::IO, O::Operator{OT,CT}) where {OT,CT}
    print(io, "Operator($(O.ops), $(O.coeffs)")
end

function Base.summary(io::IO, O::Operator)
    print(io, "Operator with $(length(O.ops)) ops and $(length(O.coeffs)) coeffs")
end

function Base.show(io::IO, ::MIME"text/plain", O::Operator{OT,CT}) where {OT,CT}
    Base.summary(io, O)
    println(io, "\n ops::Vector{$OT}:")
    for op in O.ops
        println(io, "  $op")
    end
    println(io, " coeffs: $(O.coeffs)")
end


function Base.Array{T}(O::Operator) where {T}
    A = Array{T}(O.ops[1])
    drift_offset = length(O.ops) - length(O.coeffs)
    if drift_offset == 0
        lmul!(O.coeffs[1], A)
    end
    N = length(O.ops)
    for i = 2:N
        X = Array{T}(O.ops[i])
        if i > drift_offset
            a = O.coeffs[i-drift_offset]
        else
            a = true
        end
        axpy!(a, X, A)
    end
    return A
end

Base.Array(O::Operator) = Array{eltype(O)}(O)


function Base.copy(O::Operator)
    return Operator([copy(op) for op in O.ops], copy(O.coeffs))
end


function Base.copyto!(tgt::Operator, src::Operator)
    for (i, op) in enumerate(src.ops)
        copyto!(tgt.ops[i], op)
    end
    copyto!(tgt.coeffs, src.coeffs)
end

Base.size(O::Operator) = size(O.ops[1])
Base.size(O::Operator, dim::Integer) = size(O.ops[1], dim)
Base.eltype(::Type{Operator{OT,CT}}) where {OT,CT} = promote_type(eltype(OT), CT)

function Base.getindex(O::Operator, i::Int, j::Int)
    drift_offset = length(O.ops) - length(O.coeffs)
    result = (drift_offset == 0 ? O.coeffs[1] : one(eltype(O))) * O.ops[1][i, j]
    for k = 2:length(O.ops)
        c = k <= drift_offset ? one(eltype(O)) : O.coeffs[k-drift_offset]
        result += c * O.ops[k][i, j]
    end
    return result
end

Base.length(O::Operator) = prod(size(O))

function Base.iterate(O::Operator, k = 1)
    n = length(O)
    k > n && return nothing
    n_rows = size(O, 1)
    i = (k - 1) % n_rows + 1
    j = (k - 1) ÷ n_rows + 1
    return (O[i, j], k + 1)
end


# We define `similar` to return a standard `Array` because `Operator` does not
# implement `setindex!`. So it it questionable to what extent returning an
# `Operator` would create an "uninitialized _mutable_ array". It would be
# possible to apply `similar` recursively and get something that at least we
# can `copyto!`. However, there is currently nothing in the JuliaQuantumControl
# codebase that would require this, so we're leaving that for a later time.
Base.similar(O::Operator) = Array{eltype(O)}(undef, size(O))
Base.similar(O::Operator, ::Type{S}) where {S} = Array{S}(undef, size(O))
Base.similar(O::Operator, dims::Tuple{Vararg{Int}}) = Array{eltype(O)}(undef, dims)
Base.similar(O::Operator, ::Type{S}, dims::Tuple{Vararg{Int}}) where {S} =
    Array{S}(undef, dims)


function LinearAlgebra.ishermitian(O::Operator{OT,CT}) where {OT,CT}
    return all(ishermitian(op) for op in O.ops) && all(isreal(c) for c in O.coeffs)
end


function evaluate(op::Operator, args...; kwargs...)
    return op
end


"""A static operator with a scalar pre-factor.

```julia
op = ScaledOperator(α, Ĥ)
```

represents the "lazy" product ``α Ĥ`` where ``Ĥ`` is an operator (typically an
[`Operator`](@ref) instance) and ``α`` is a scalar.
"""
struct ScaledOperator{CT<:Number,OT}
    coeff::CT
    operator::OT

    function ScaledOperator(coeff::CT, operator::OT) where {CT,OT}
        if coeff == 1.0
            return operator
        else
            return new{CT,OT}(coeff, operator)
        end
    end
end

function Base.summary(io::IO, O::ScaledOperator)
    print(io, "ScaledOperator with coeff=$(O.coeff)")
end

function Base.show(io::IO, O::ScaledOperator{CT,OT}) where {CT,OT}
    print(io, "ScaledOperator($(O.coeff), $(O.operator))")
end

function Base.show(io::IO, ::MIME"text/plain", O::ScaledOperator{CT,OT}) where {CT,OT}
    Base.summary(io, O)
    println(io, "\n operator.ops::$(typeof(O.operator.ops)):")
    for op in O.operator.ops
        println(io, "  $op")
    end
    println(io, " operator.coeffs: $(O.operator.coeffs)")
end


function Base.Array{T}(O::ScaledOperator{CT,Operator{OOT,OCT}}) where {T,CT<:Number,OOT,OCT}
    A = Array{T}(O.operator.ops[1])
    drift_offset = length(O.operator.ops) - length(O.operator.coeffs)
    a = O.coeff
    (drift_offset == 0) && (a *= O.operator.coeffs[1])
    lmul!(a, A)
    N = length(O.operator.ops)
    for i = 2:N
        X = Array{T}(O.operator.ops[i])
        a = O.coeff
        (i > drift_offset) && (a *= O.operator.coeffs[i-drift_offset])
        axpy!(a, X, A)
    end
    return A
end

# fallback (less efficient, but doesn't assume Operator-OT)
Base.Array{T}(O::ScaledOperator{CT,OT}) where {T,CT,OT} = O.coeff * Array{T}(O.operator)
Base.Array(O::ScaledOperator) = Array{eltype(O)}(O)


Base.size(O::ScaledOperator) = size(O.operator)
Base.size(O::ScaledOperator, dim::Integer) = size(O.operator, dim)
Base.eltype(::Type{ScaledOperator{CT,OT}}) where {CT,OT} = promote_type(CT, eltype(OT))

function Base.getindex(O::ScaledOperator, i::Int, j::Int)
    return O.coeff * O.operator[i, j]
end

Base.length(O::ScaledOperator) = prod(size(O))

function Base.iterate(O::ScaledOperator, k = 1)
    n = length(O)
    k > n && return nothing
    n_rows = size(O, 1)
    i = (k - 1) % n_rows + 1
    j = (k - 1) ÷ n_rows + 1
    return (O[i, j], k + 1)
end

Base.similar(O::ScaledOperator) = Array{eltype(O)}(undef, size(O))
Base.similar(O::ScaledOperator, ::Type{S}) where {S} = Array{S}(undef, size(O))
Base.similar(O::ScaledOperator, dims::Tuple{Vararg{Int}}) = Array{eltype(O)}(undef, dims)
Base.similar(O::ScaledOperator, ::Type{S}, dims::Tuple{Vararg{Int}}) where {S} =
    Array{S}(undef, dims)

LinearAlgebra.ishermitian(O::ScaledOperator) = (isreal(O.coeff) && ishermitian(O.operator))


function evaluate(op::ScaledOperator, args...; kwargs...)
    return op
end


"""Initialize a (usually time-dependent) Hamiltonian.

The most common usage is, e.g.,

```jldoctest hamiltonian
using QuantumPropagators

H₀ = ComplexF64[0 0; 0 1];
H₁ = ComplexF64[0 1; 1 0];
ϵ₁(t) = 1.0;

hamiltonian(H₀, (H₀, ϵ₁))

# output

Generator with 2 ops and 1 amplitudes
 ops::Vector{Matrix{ComplexF64}}:
  ComplexF64[0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 1.0 + 0.0im]
  ComplexF64[0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 1.0 + 0.0im]
 amplitudes::Vector{typeof(ϵ₁)}:
  ϵ₁
```

In general,

```julia
H = hamiltonian(terms...; check=true)
```

constructs a Hamiltonian based on the given `terms`. Each term must be an
operator or a tuple `(op, ampl)` of an operator and a control amplitude. Single
operators are considered "drift" terms.

In most cases, each control amplitude will simply be a control function or
vector of pulse values. In general, `ampl` can be an arbitrary object that
depends on one or more controls, which must be obtainable via
[`get_controls(ampl)`](@ref get_controls). See
`QuantumPropagators.Interfaces.check_amplitude` for the required
interface.

The `hamiltonian` function will generally return a [`Generator`](@ref)
instance. However, if none of the given terms are time-dependent, it may return
a static operator (e.g., an `AbstractMatrix` or [`Operator`](@ref)):

```jldoctest hamiltonian
hamiltonian(H₀)
# output
2×2 Matrix{ComplexF64}:
 0.0+0.0im  0.0+0.0im
 0.0+0.0im  1.0+0.0im
```

```jldoctest hamiltonian
hamiltonian(H₀, (H₁, 2.0))
# output
Operator with 2 ops and 1 coeffs
 ops::Vector{Matrix{ComplexF64}}:
  ComplexF64[0.0 + 0.0im 0.0 + 0.0im; 0.0 + 0.0im 1.0 + 0.0im]
  ComplexF64[0.0 + 0.0im 1.0 + 0.0im; 1.0 + 0.0im 0.0 + 0.0im]
 coeffs: [2.0]
```

The `hamiltonian` function may generate warnings if the `terms` are of an
unexpected type or structure.  These can be suppressed with `check=false`.
"""
hamiltonian(terms...; check = true) = _make_generator(terms...; check)

function _make_generator(terms...; check = false)
    ops = Any[]
    drift = Any[]
    amplitudes = Any[]
    if check
        if (length(terms) == 1) && (terms[1] isa Union{Tuple,Vector})
            @warn("Generator terms may not have been properly expanded")
        end
    end
    for term in terms
        if term isa Union{Tuple,Vector}
            length(term) == 2 || error("time-dependent term must be 2-tuple")
            op, ampl = term
            if check
                StdAmplType = Union{Function,Vector}
                if (op isa StdAmplType) && !(ampl isa StdAmplType)
                    @warn("It looks like (op, ampl) in term are reversed")
                end
            end
            i = findfirst(a -> a == ampl, amplitudes)
            if isnothing(i)
                push!(ops, op)
                push!(amplitudes, ampl)
            else
                if check
                    try
                        ops[i] = ops[i] + op
                    catch exc
                        @error(
                            "Collected operators are of a disparate type: $(typeof(ops[i])), $(typeof(op)): $exc"
                        )
                        rethrow()
                    end
                else
                    ops[i] = ops[i] + op
                end
            end
        else
            op = term
            if length(drift) == 0
                push!(drift, op)
            else
                if check
                    try
                        drift[1] = drift[1] + op
                    catch exc
                        @error(
                            "Collected drift operators are of a disparate type: $(typeof(drift[1])), $(typeof(op)): $exc"
                        )
                        rethrow()
                    end
                else
                    drift[1] = drift[1] + op
                end
            end
        end
    end
    ops = [drift..., ops...]  # narrow eltype
    OT = eltype(ops)
    amplitudes = [amplitudes...]  # narrow eltype
    AT = eltype(amplitudes)
    if length(amplitudes) == 0
        (length(drift) > 0) || error("Generator has no terms")
        return drift[1]
    else
        if check
            if !(isconcretetype(OT))
                @warn("Collected operators are not of a concrete type: $OT")
            end
            if AT ≡ Any
                @warn("Collected amplitudes are of disparate types")
            end
        end
        if (AT <: Number)
            return Operator(ops, amplitudes)
        else
            return Generator(ops, amplitudes)
        end
    end
end



function ham_to_superop(H::AbstractSparseMatrix; convention)
    # See https://arxiv.org/abs/1312.0111, Appendix B.2
    ⊗(A, B) = kron(A, B)
    𝟙 = SparseMatrixCSC{ComplexF64,Int64}(sparse(I, size(H)[1], size(H)[2]))
    H_T = sparse(transpose(H))
    L = sparse(𝟙 ⊗ H - H_T ⊗ 𝟙)
    if convention == :TDSE
        return L
    elseif convention == :LvN
        return 1im * L
    else
        throw(ArgumentError("convention must be :TDSE or :LvN"))
    end
end

function ham_to_superop(H::AbstractMatrix; convention)
    return ham_to_superop(sparse(H); convention = convention)
end


function lindblad_to_superop(A::AbstractSparseMatrix; convention)
    # See https://arxiv.org/abs/1312.0111, Appendix B.2
    ⊗(A, B) = kron(A, B)
    A⁺ = sparse(A')
    A⁺ᵀ = sparse(transpose(A⁺))
    A⁺_A = sparse(A⁺ * A)
    A⁺_A_ᵀ = sparse(transpose(A⁺_A))
    𝟙 = SparseMatrixCSC{ComplexF64,Int64}(sparse(I, size(A)[1], size(A)[2]))
    D = sparse(A⁺ᵀ ⊗ A - (𝟙 ⊗ A⁺_A) / 2 - (A⁺_A_ᵀ ⊗ 𝟙) / 2)
    if convention == :TDSE
        return 1im * D
    elseif convention == :LvN
        return D
    else
        throw(ArgumentError("convention must be :TDSE or :LvN"))
    end
end

function lindblad_to_superop(A::AbstractMatrix; convention)
    return lindblad_to_superop(sparse(A); convention = convention)
end


function dissipator(c_ops; convention)
    N = size(c_ops[1])[1]
    @assert N == size(c_ops[1])[2]
    D = spzeros(ComplexF64, N^2, N^2)
    for A in c_ops
        D += lindblad_to_superop(A; convention = convention)
    end
    return (D,)
end


nhilbert(H::AbstractMatrix) = size(H)[1]
nhilbert(H::Tuple{HT,ET}) where {HT<:AbstractMatrix,ET} = size(H[1])[1]


@doc raw"""Construct a Liouvillian [`Generator`](@ref).

```julia
ℒ = liouvillian(Ĥ, c_ops=(); convention=:LvN, check=true)
```

calculates the sparse Liouvillian super-operator `ℒ` from the Hamiltonian `Ĥ`
and a list `c_ops` of Lindblad operators.

With `convention=:LvN`, applying the resulting `ℒ` to a vectorized density
matrix `ρ⃗` calculates ``\frac{d}{dt} \vec{\rho}(t) = ℒ \vec{\rho}(t)``
equivalent to the Liouville-von-Neumann equation for the density matrix ``ρ̂``,

```math
\frac{d}{dt} ρ̂(t)
= -i [Ĥ, ρ̂(t)] + \sum_k\left(
    Â_k ρ̂ Â_k^\dagger
    - \frac{1}{2} A_k^\dagger Â_k ρ̂
    - \frac{1}{2} ρ̂ Â_k^\dagger Â_k
  \right)\,,
```

where the Lindblad operators ``Â_k`` are the elements of `c_ops`.

The Hamiltonian ``Ĥ`` will generally be time-dependent. For example, it may be
a [`Generator`](@ref) as returned by [`hamiltonian`](@ref). For example, for a
Hamiltonian with the terms `(Ĥ₀, (Ĥ₁, ϵ₁), (Ĥ₂, ϵ₂))`, where `Ĥ₀`, `Ĥ₁`, `Ĥ₂`
are matrices, and `ϵ₁` and `ϵ₂` are functions of time, the resulting `ℒ` will
be a [`Generator`](@ref) corresponding to terms `(ℒ₀, (ℒ₁, ϵ₁), (ℒ₂, ϵ₂))`,
where the initial terms is the superoperator `ℒ₀` for the static component of
the Liouvillian, i.e., the commutator with the drift Hamiltonian `Ĥ₀`, plus the
dissipator (sum over ``k``), as a sparse matrix. Time-dependent Lindblad
operators are not currently supported. The remaining elements are tuples `(ℒ₁,
ϵ₁)` and `(ℒ₂, ϵ₂)` corresponding to the commutators with the two control
Hamiltonians, where `ℒ₁` and `ℒ₂` again are sparse matrices.

If ``Ĥ`` is not time-dependent, the resulting `ℒ` will likewise be a static
operator. Passing `H=nothing` with non-empty `c_ops` initializes a pure
dissipator.

With `convention=:TDSE`, the Liouvillian will be constructed for the equation
of motion ``i \hbar \frac{d}{dt} \vec{\rho}(t) = ℒ \vec{\rho}(t)`` to match
exactly the form of the time-dependent Schrödinger equation. While this
notation is not standard in the literature of open quantum systems, it has the
benefit that the resulting `ℒ` can be used in a numerical propagator for a
(non-Hermitian) Schrödinger equation without any change. Thus, for numerical
applications, `convention=:TDSE` is generally preferred. The returned `ℒ`
between the two conventions differs only by a factor of ``i``, since we
generally assume ``\hbar=1``.

The `convention` keyword argument is mandatory, to force a conscious choice.

See [Goerz et. al. "Optimal control theory for a unitary operation under
dissipative evolution", arXiv 1312.0111v2, Appendix
B.2](https://arxiv.org/abs/1312.0111v2) for the explicit construction of the
Liouvillian superoperator as a sparse matrix.

Passing `check=false`, suppresses warnings and errors about unexpected types or
the structure of the arguments, cf. [`hamiltonian`](@ref).
"""
function liouvillian(H::Tuple, c_ops = (); kwargs...)
    check = get(kwargs, :check, true)
    return liouvillian(_make_generator(H...; check), c_ops; kwargs...)
end

function liouvillian(H::Generator, c_ops = (); convention, check = true)
    terms = []
    if length(c_ops) > 0
        append!(terms, dissipator(c_ops; convention))
    end
    drift_offset = length(H.ops) - length(H.amplitudes)
    for (i, op) in enumerate(H.ops)
        if i <= drift_offset
            term = ham_to_superop(op; convention)
            push!(terms, term)
        else
            ampl = H.amplitudes[i-drift_offset]
            term = (ham_to_superop(op; convention), ampl)
            push!(terms, term)
        end
    end
    return _make_generator(terms...; check)
end

function liouvillian(H::AbstractMatrix, c_ops = (); convention, check = true)
    L0 = ham_to_superop(H; convention)
    terms = Any[L0,]
    if length(c_ops) > 0
        append!(terms, dissipator(c_ops; convention))
    end
    return _make_generator(terms...; check)
end

function liouvillian(H::Nothing, c_ops = (); convention, check = true)
    if length(c_ops) > 0
        terms = dissipator(c_ops; convention = convention)
        return _make_generator(terms...; check)
    else
        error("Empty Liouvillian, must give at least one of `H` or `c_ops`")
    end
end


function LinearAlgebra.mul!(C, A::Operator, B, α, β)
    drift_offset = length(A.ops) - length(A.coeffs)
    c = α
    (drift_offset == 0) && (c *= A.coeffs[1])
    mul!(C, A.ops[1], B, c, β)
    for i = 2:length(A.ops)
        c = α
        (i > drift_offset) && (c *= A.coeffs[i-drift_offset])
        mul!(C, A.ops[i], B, c, true)
    end
    return C
end


function LinearAlgebra.dot(x, A::Operator, y)
    drift_offset = length(A.ops) - length(A.coeffs)
    result::ComplexF64 = 0
    for i = 1:length(A.ops)
        if i > drift_offset
            c = A.coeffs[i-drift_offset]
            result += c * dot(x, A.ops[i], y)
        else
            result += dot(x, A.ops[i], y)
        end
    end
    return result
end


function Base.:*(α::Number, O::Operator)
    return ScaledOperator(α, O)
end

Base.:*(O::Operator, α::Number) = α * O


# Not-in-place |Φ⟩ = c * Ô * |Ψ⟩
@inline function _op_mul_psi(O::Operator, Ψ, c)
    drift_offset = length(O.ops) - length(O.coeffs)
    α = c
    (drift_offset == 0) && (α *= O.coeffs[1])
    Φ = α * O.ops[1] * Ψ
    for i = 2:length(O.ops)
        α = c
        (i > drift_offset) && (α *= O.coeffs[i-drift_offset])
        Φ = Φ + α * O.ops[i] * Ψ
    end
    return Φ
end

Base.:*(O::Operator, Ψ) = _op_mul_psi(O, Ψ, 1.0)


Base.convert(::Type{MT}, O::Operator) where {MT<:Matrix} = convert(MT, Array(O))
Base.convert(::Type{MT}, O::ScaledOperator) where {MT<:Matrix} = convert(MT, Array(O))


function Base.:*(α::Number, O::ScaledOperator)
    return ScaledOperator(α * O.coeff, O.operator)
end

Base.:*(O::ScaledOperator, α::Number) = α * O


Base.:*(O::ScaledOperator, Ψ) = _op_mul_psi(O.operator, Ψ, O.coeff)


function LinearAlgebra.mul!(C, A::ScaledOperator, B, α, β)
    return mul!(C, A.operator, B, A.coeff * α, β)
end


function LinearAlgebra.dot(x, A::ScaledOperator, y)
    return A.coeff * dot(x, A.operator, y)
end


function get_controls(generator::Generator)
    controls = []
    slots_dict = IdDict()  # utilized as Set of controls we've seen
    for (i, ampl) in enumerate(generator.amplitudes)
        for control in get_controls(ampl)
            if control in keys(slots_dict)
                # We've seen this control before, so we just record the slot
                # where it is referenced
                push!(slots_dict[control], i)
            else
                push!(controls, control)
                slots_dict[control] = [i]
            end
        end
    end
    return Tuple(controls)
end


function get_controls(generator::Tuple)
    G = _make_generator(generator...)
    return get_controls(G)
end


get_controls(operator::Operator) = Tuple([])
get_controls(operator::ScaledOperator) = Tuple([])


function evaluate(generator::Generator, args...; vals_dict = IdDict())
    coeffs = []
    for (i, ampl) in enumerate(generator.amplitudes)
        coeff = evaluate(ampl, args...; vals_dict)
        if coeff isa Number
            push!(coeffs, coeff)
        else
            error(
                "In `evaluate($generator, $args, vals_dict=$vals_dict)`, the amplitude $i evaluates to $(typeof(coeff)), not a number"
            )
        end
    end
    coeffs = [coeffs...]  # narrow eltype
    return Operator(generator.ops, coeffs)
end


function evaluate!(op::Operator, generator::Generator, args...; vals_dict = IdDict())
    @assert length(op.ops) == length(generator.ops)
    @assert all(O ≡ P for (O, P) in zip(op.ops, generator.ops))
    for (i, ampl) in enumerate(generator.amplitudes)
        coeff = evaluate(ampl, args...; vals_dict)
        @assert coeff isa Number
        op.coeffs[i] = coeff
    end
    return op
end


function substitute(generator::Generator, replacements)
    if generator ∈ keys(replacements)
        return replacements[generator]
    end
    ops = [substitute(op, replacements) for op in generator.ops]
    amplitudes = [substitute(ampl, replacements) for ampl in generator.amplitudes]
    return Generator(ops, amplitudes)
end

function substitute(operator::Operator, replacements)
    if operator ∈ keys(replacements)
        return replacements[operator]
    end
    ops = [substitute(op, replacements) for op in operator.ops]
    return Operator(ops, operator.coeffs)
end


end
