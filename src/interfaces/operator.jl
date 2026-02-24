using Test

using LinearAlgebra
import ArrayInterface
using ..Controls: get_controls, evaluate


"""Check that `op` is a valid operator that can be applied to `state`.

```julia
@test check_operator(op; state, tlist=[0.0, 1.0],
                     for_expval=true, atol=1e-14, quiet=false)
```

verifies the given `op` relative to `state`. The `state` must pass
[`check_state`](@ref).

An "operator" is any object that [`evaluate`](@ref) returns when evaluating a
time-dependent dynamic generator. The specific requirements for `op` are:

* `size(op)` must be defined and return a tuple of integers
* `size(op, dim)` must be defined for each dimension and be consistent with
  `size(op)`
* `op` must not be time-dependent: [`evaluate(op, tlist, 1) ≡ op`](@ref evaluate)
* `op` must not contain any controls: [`length(get_controls(op)) == 0`](@ref get_controls)
* `op * state` must be defined
* The [`QuantumPropagators.Interfaces.supports_inplace`](@ref) method must be
  defined for `op`. If it returns `true`, it must be possible to evaluate a
  generator in-place into the existing `op`. See
  [`QuantumPropagators.Interfaces.check_generator`](@ref).

If [`QuantumPropagators.Interfaces.supports_inplace(state)`](@ref
QuantumPropagators.Interfaces.supports_inplace):

* The 3-argument `LinearAlgebra.mul!` must apply `op` to the given `state`
* The 5-argument `LinearAlgebra.mul!` must apply `op` to the given `state`
* `LinearAlgebra.mul!` must match `*`, if applicable
* `LinearAlgebra.mul!` must return the resulting state

If `for_expval` (typically required for optimal control):

* `LinearAlgebra.dot(state, op, state)` must return return a number
* `dot(state, op, state)` must match `dot(state, op * state)`, if applicable

If [`QuantumPropagators.Interfaces.supports_matrix_interface(op)`](@ref
QuantumPropagators.Interfaces.supports_matrix_interface) is `true`, the
operator must implement the
[Abstract Array interface](https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-array)
for two-dimensional arrays:

* `eltype(op)` must be defined and return a numeric type
* `getindex(op, i, j)` must be defined and return elements matching `eltype`
* `length(op)` must equal `prod(size(op))`
* `iterate(op)` must be defined
* `similar(op)` must be defined and return a mutable array with the same shape
  and element type. "Mutability" is determined by
  [`ArrayInterface.ismutable`](@extref).
* `similar(op, ::Type{S})` must return a mutable array with the same shape and
  element type `S`
* `similar(op, dims::Dims)` must return a mutable array with the same element
  type and the given dimensions
* `similar(op, ::Type{S}, dims::Dims)` must return a mutable array with the
  given element type and dimensions

The read-write method `setindex!` is not a requirement for
`supports_matrix_interface`.

The function returns `true` for a valid operator and `false` for an invalid
operator. Unless `quiet=true`, it will log an error to indicate which of the
conditions failed.
"""
function check_operator(
    op;
    state,
    tlist = [0.0, 1.0],
    for_expval = true,
    atol = 1e-14,
    quiet = false,
    _message_prefix = ""  # for recursive calling
)

    ≈(a, b) = isapprox(a, b; atol)

    px = _message_prefix
    success = true

    Ψ = state

    @assert check_state(state; atol, quiet = true)
    @assert tlist isa Vector{Float64}

    try
        supports_inplace(op)
    catch exc
        quiet || @error(
            "$(px)The `QuantumPropagators.Interfaces.supports_inplace` method must be defined for type `$(typeof(op))`.",
            exception = (exc, catch_abbreviated_backtrace())
        )
        success = false
    end

    try
        s = size(op)
        if !(s isa Tuple)
            quiet || @error "$(px)`size(op)` must return a tuple, not $(typeof(s))"
            success = false
        elseif !all(d isa Integer for d in s)
            quiet || @error "$(px)`size(op)` must return a tuple of integers, not $s"
            success = false
        else
            for (dim, n) in enumerate(s)
                try
                    n_dim = size(op, dim)
                    if !(n_dim isa Integer)
                        quiet ||
                            @error "$(px)`size(op, $dim)` must return an integer, not $(typeof(n_dim))"
                        success = false
                    elseif n_dim != n
                        quiet ||
                            @error "$(px)`size(op, $dim)` must be consistent with `size(op)`: $n_dim ≠ $n"
                        success = false
                    end
                catch exc
                    quiet || @error(
                        "$(px)`size(op, $dim)` must be defined.",
                        exception = (exc, catch_abbreviated_backtrace())
                    )
                    success = false
                end
            end
        end
    catch exc
        quiet || @error(
            "$(px)`size(op)` must be defined.",
            exception = (exc, catch_abbreviated_backtrace())
        )
        success = false
    end

    try
        H = evaluate(op, tlist, 1)
        if H ≢ op
            quiet || @error "$(px)`evaluate(op, tlist, 1) must return operator ≡ op"
            success = false
        end
    catch exc
        quiet || @error(
            "$(px)op must not be time-dependent.",
            exception = (exc, catch_abbreviated_backtrace())
        )
        success = false
    end

    try
        controls = get_controls(op)
        if length(controls) ≠ 0
            quiet ||
                @error "$(px)get_controls(op) must return an empty tuple, not $controls"
            success = false
        end
    catch exc
        quiet || @error(
            "$(px)op must not contain any controls.",
            exception = (exc, catch_abbreviated_backtrace())
        )
        success = false
    end


    try
        ϕ = op * Ψ
        if !(ϕ isa typeof(Ψ))
            quiet ||
                @error "$(px)`op * state` must return an object of the same type as `state`, not $typeof(ϕ)"
            success = false
        end
    catch exc
        quiet || @error(
            "$(px)`op * state` must be defined.",
            exception = (exc, catch_abbreviated_backtrace())
        )
        success = false
    end

    if supports_inplace(state)

        try
            ϕ = similar(Ψ)
            if mul!(ϕ, op, Ψ) ≢ ϕ
                quiet || @error "$(px)`mul!(ϕ, op, state)` must return the resulting ϕ"
                success = false
            end
            Δ = norm(ϕ - op * Ψ)
            if Δ > atol
                quiet || @error "$(px)`mul!(ϕ, op, state)` must match `op * state`" Δ atol
                success = false
            end
        catch exc
            quiet || @error(
                "$(px)The 3-argument `mul!` must apply `op` to the given `state`.",
                exception = (exc, catch_abbreviated_backtrace())
            )
            success = false
        end

        try
            ϕ = similar(Ψ)
            ϕ0 = similar(Ψ)
            copyto!(ϕ, Ψ)
            copyto!(ϕ0, Ψ)
            if mul!(ϕ, op, Ψ, 0.5, 0.5) ≢ ϕ
                quiet ||
                    @error "$(px)`mul!(ϕ, op, state, α, β)` must return the resulting ϕ"
                success = false
            end
            Δ = norm(ϕ - (0.5 * ϕ0 + 0.5 * (op * Ψ)))
            if Δ > atol
                quiet ||
                    @error "$(px)`mul!(ϕ, op, state, α, β)` must match β*ϕ + α*op*state" Δ atol
                success = false
            end
        catch exc
            quiet || @error(
                "$(px)The 5-argument `mul!` must apply `op` to the given `state`.",
                exception = (exc, catch_abbreviated_backtrace())
            )
            success = false
        end

    end

    if for_expval

        try
            val = dot(Ψ, op, Ψ)
            if !(val isa Number)
                quiet ||
                    @error "$(px)`dot(state, op, state)` must return a number, not $typeof(val)"
                success = false
            end
            Δ = abs(dot(Ψ, op, Ψ) - dot(Ψ, op * Ψ))
            if Δ > atol
                quiet ||
                    @error "$(px)`dot(state, op, state)` must match `dot(state, op * state)`" Δ atol
                success = false
            end
            if supports_inplace(state)
                ϕ = similar(Ψ)
                mul!(ϕ, op, Ψ)
                Δ = abs(dot(Ψ, op, Ψ) - dot(Ψ, ϕ))
                if Δ > atol
                    quiet ||
                        @error "$(px)`dot(state, op, state)` must match `dot(state, op * state)`" Δ atol
                    success = false
                end
            end
        catch exc
            quiet || @error(
                "$(px)`dot(state, op, state)` must return return a number.",
                exception = (exc, catch_abbreviated_backtrace())
            )
            success = false
        end

    end

    if supports_matrix_interface(op)

        try
            T = eltype(op)
            if !(T isa Type && T <: Number)
                quiet || @error "$(px)`eltype(op)` must return a numeric type, not $T"
                success = false
            end
        catch exc
            quiet || @error(
                "$(px)`eltype(op)` must be defined.",
                exception = (exc, catch_abbreviated_backtrace())
            )
            success = false
        end

        try
            s = size(op)
            if length(s) == 2 && all(d -> d > 0, s)
                val = op[1, 1]
                T = eltype(op)
                if T isa Type && T <: Number && !(val isa T)
                    quiet ||
                        @error "$(px)`op[1, 1]` must return a value of type `eltype(op)=$T`, not $(typeof(val))"
                    success = false
                end
            end
        catch exc
            quiet || @error(
                "$(px)`getindex(op, i, j)` must be defined.",
                exception = (exc, catch_abbreviated_backtrace())
            )
            success = false
        end

        try
            l = length(op)
            s = size(op)
            if l != prod(s)
                quiet ||
                    @error "$(px)`length(op)` must equal `prod(size(op))`: $l ≠ $(prod(s))"
                success = false
            end
        catch exc
            quiet || @error(
                "$(px)`length(op)` must be defined.",
                exception = (exc, catch_abbreviated_backtrace())
            )
            success = false
        end

        try
            itr = iterate(op)
            s = size(op)
            if isnothing(itr) && prod(s) > 0
                quiet ||
                    @error "$(px)`iterate(op)` must not return `nothing` for a non-empty operator"
                success = false
            end
        catch exc
            quiet || @error(
                "$(px)`iterate(op)` must be defined.",
                exception = (exc, catch_abbreviated_backtrace())
            )
            success = false
        end

        try
            op2 = similar(op)
            if !ArrayInterface.ismutable(op2)
                quiet ||
                    @error "$(px)`similar(op)` must return a mutable array (`ArrayInterface.ismutable` must be `true`), got $(typeof(op2))"
                success = false
            end
            if size(op2) != size(op)
                quiet ||
                    @error "$(px)`similar(op)` must return an array with the same shape: size $(size(op2)) ≠ $(size(op))"
                success = false
            end
            if eltype(op2) != eltype(op)
                quiet ||
                    @error "$(px)`similar(op)` must return an array with the same element type: $(eltype(op2)) ≠ $(eltype(op))"
                success = false
            end
        catch exc
            quiet || @error(
                "$(px)`similar(op)` must be defined.",
                exception = (exc, catch_abbreviated_backtrace())
            )
            success = false
        end

        try
            S = (eltype(op) == ComplexF64) ? ComplexF32 : ComplexF64
            op2 = similar(op, S)
            if !ArrayInterface.ismutable(op2)
                quiet ||
                    @error "$(px)`similar(op, $S)` must return a mutable array (`ArrayInterface.ismutable` must be `true`), got $(typeof(op2))"
                success = false
            end
            if size(op2) != size(op)
                quiet ||
                    @error "$(px)`similar(op, $S)` must return an array with the same shape: size $(size(op2)) ≠ $(size(op))"
                success = false
            end
            if eltype(op2) != S
                quiet ||
                    @error "$(px)`similar(op, $S)` must return an array with element type $S, got $(eltype(op2))"
                success = false
            end
        catch exc
            quiet || @error(
                "$(px)`similar(op, ::Type{S})` must be defined.",
                exception = (exc, catch_abbreviated_backtrace())
            )
            success = false
        end

        try
            dims = size(op)
            op2 = similar(op, dims)
            if !ArrayInterface.ismutable(op2)
                quiet ||
                    @error "$(px)`similar(op, dims)` must return a mutable array (`ArrayInterface.ismutable` must be `true`), got $(typeof(op2))"
                success = false
            end
            if size(op2) != dims
                quiet ||
                    @error "$(px)`similar(op, dims)` must return an array with size $dims, got $(size(op2))"
                success = false
            end
            if eltype(op2) != eltype(op)
                quiet ||
                    @error "$(px)`similar(op, dims)` must return an array with the same element type: $(eltype(op2)) ≠ $(eltype(op))"
                success = false
            end
        catch exc
            quiet || @error(
                "$(px)`similar(op, dims::Dims)` must be defined.",
                exception = (exc, catch_abbreviated_backtrace())
            )
            success = false
        end

        try
            S = (eltype(op) == ComplexF64) ? ComplexF32 : ComplexF64
            dims = size(op)
            op2 = similar(op, S, dims)
            if !ArrayInterface.ismutable(op2)
                quiet ||
                    @error "$(px)`similar(op, $S, dims)` must return a mutable array (`ArrayInterface.ismutable` must be `true`), got $(typeof(op2))"
                success = false
            end
            if size(op2) != dims
                quiet ||
                    @error "$(px)`similar(op, $S, dims)` must return an array with size $dims, got $(size(op2))"
                success = false
            end
            if eltype(op2) != S
                quiet ||
                    @error "$(px)`similar(op, $S, dims)` must return an array with element type $S, got $(eltype(op2))"
                success = false
            end
        catch exc
            quiet || @error(
                "$(px)`similar(op, ::Type{S}, dims::Dims)` must be defined.",
                exception = (exc, catch_abbreviated_backtrace())
            )
            success = false
        end

    end

    return success

end
