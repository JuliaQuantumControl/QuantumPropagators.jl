using Test

using LinearAlgebra


"""Check that `state` is a valid element of a Hilbert space.

```julia
@test check_state(state; normalized=false, atol=1e-15, quiet=false)
```

verifies the following requirements:

* The inner product (`LinearAlgebra.dot`) of two states must return a Complex
  number.
* The `LinearAlgebra.norm` of `state` must be defined via the inner product.
  This is the *definition* of a Hilbert space, a.k.a a "complete inner product
  space" or more precisely a "Banach space (normed vector space) where the
  norm is induced by an inner product".
* The [`QuantumPropagators.Interfaces.supports_inplace](@ref) method must be
  defined for `state`

Any `state` must support the following not-in-place operations:

* `state + state` and `state - state` must be defined
* `copy(state)` must be defined and return an object of the same type as
  `state`
* `c * state` for a scalar `c` must be defined
* `norm(state + state)` must fulfill the triangle inequality
* `zero(state)` must be defined and produce a state with norm 0
* `0.0 * state` must produce a state with norm 0
* `copy(state) - state` must have norm 0
* `norm(state)` must have absolute homogeneity: `norm(s * state) = s *
  norm(state)`

If `supports_inplace(state)` is `true`, the `state` must also support the
following:

* `similar(state)` must be defined and return a valid state of the same type a
  `state`
* `copyto!(other, state)` must be defined
* `fill!(state, c)` must be defined
* `LinearAlgebra.lmul!(c, state)` for a scalar `c` must be defined
* `LinearAlgebra.axpy!(c, state, other)` must be defined
* `norm(state)` must fulfill the same general mathematical norm properties as
  for the non-in-place norm.

If `normalized` (not required by default):

* `LinearAlgebra.norm(state)` must be 1

If [`QuantumPropagators.Interfaces.supports_vector_interface(state)`](@ref
QuantumPropagators.Interfaces.supports_vector_interface) is `true`, the state
must implement the
[Abstract Array interface](https://docs.julialang.org/en/v1/manual/interfaces/#man-interface-array)
for one-dimensional arrays:

* `eltype(state)` must be defined and return a numeric type
* `getindex(state, i)` must be defined and return elements matching `eltype`
* `length(state)` must equal `prod(size(state))`
* `iterate(state)` must be defined
* `similar(state)` must be defined and return a mutable vector with the same
  length and element type.
* `similar(state, ::Type{S})` must return a mutable vector with the same length
  and element type `S`
* `similar(state, dims::Dims)` must return a mutable array with the same
  element type and the given dimensions
* `similar(state, ::Type{S}, dims::Dims)` must return a mutable array with the
  given element type and dimensions

If additionally [`QuantumPropagators.Interfaces.supports_inplace(state)`](@ref
QuantumPropagators.Interfaces.supports_inplace) is `true` (read-write vector):

* `setindex!(state, v, i)` must be defined

Note: When `supports_inplace(state)` is `true`, `similar(state)` is already
required to return the *same type* as `state` (see above). The vector interface
checks are weaker and complementary.

It is strongly recommended to always support immutable operations (also for
mutable states)

The function returns `true` for a valid state and `false` for an invalid state.
Unless `quiet=true`, it will log an error to indicate which of the conditions
failed.
"""
function check_state(
    state;
    normalized = false,
    atol = 1e-14,
    quiet = false,
    _check_similar = true,  # to avoid infinite recursion
    _message_prefix = ""  # for recursive calling
)

    ≈(a, b) = isapprox(a, b; atol)
    px = _message_prefix

    success = true
    inplace = false

    try
        inplace = supports_inplace(state)
    catch exc
        quiet || @error(
            "$(px)The `QuantumPropagators.Interfaces.supports_inplace` method must be defined for type `$(typeof(state))`.",
            exception = (exc, catch_abbreviated_backtrace())
        )
        success = false
    end

    try
        c = state ⋅ state
        if !(c isa Complex)
            quiet ||
                @error "$(px)`state ⋅ state` must return a Complex number type, not $(typeof(c))"
            success = false
        end
    catch exc
        quiet || @error(
            "$(px)the inner product of two states must be a complex number.",
            exception = (exc, catch_abbreviated_backtrace())
        )
        success = false
    end

    try
        r1 = norm(state)
        r2 = sqrt(real(state ⋅ state))
        Δ = abs(r1 - r2)
        if Δ > atol
            quiet || @error "$(px)`norm(state)=$r1)` must match `√(state⋅state)=$r2`" Δ atol
            success = false
        end
    catch exc
        quiet || @error(
            "$(px)the norm of a state must be defined via the inner product.",
            exception = (exc, catch_abbreviated_backtrace())
        )
        success = false
    end

    # not-in-place interface:

    try
        Δ = norm(state - state)
        if Δ > atol
            quiet || @error "`$(px)state - state` must have norm 0" Δ atol
            success = false
        end
        η = norm(state + state)
        if η > (2 * norm(state) + atol)
            quiet ||
                @error "`$(px)norm(state + state)` must fulfill the triangle inequality"
            success = false
        end
    catch exc
        quiet || @error(
            "$(px)`state + state` and `state - state` must be defined.",
            exception = (exc, catch_abbreviated_backtrace())
        )
        success = false
    end

    try
        ϕ = copy(state)
        if typeof(ϕ) != typeof(state)
            quiet ||
                @error "$(px)`copy(state)::$(typeof(ϕ))` must have the same type as `state::$(typeof(state))`"
            success = false
        end
        Δ = norm(ϕ - state)
        if Δ > atol
            quiet || @error "$(px)`copy(state) - state` must have norm 0" Δ atol
            success = false
        end
    catch exc
        quiet || @error(
            "$(px)copy(state) must be defined.",
            exception = (exc, catch_abbreviated_backtrace())
        )
        success = false
    end

    try
        ϕ = 0.5 * state
        Δ = abs(norm(ϕ) - 0.5 * norm(state))
        if Δ > atol
            quiet ||
                @error "$(px)`norm(state)` must have absolute homogeneity: `norm(s * state) = s * norm(state)`" Δ atol
            success = false
        end
    catch exc
        quiet || @error(
            "$(px)`c * state` for a scalar `c` must be defined.",
            exception = (exc, catch_abbreviated_backtrace())
        )
        success = false
    end

    try
        state_zero = zero(state)
        Δ = norm(state_zero)
        if Δ > atol
            quiet || @error "$(px)`zero(state)` must produce a state with norm 0." Δ atol
        end
    catch exc
        quiet || @error(
            "$(px)`zero(state)` must be defined and produce a state with norm 0.",
            exception = (exc, catch_abbreviated_backtrace())
        )
        success = false
    end

    try
        Δ = norm(0.0 * state)
        if Δ > atol
            quiet || @error "$(px)`0.0 * state` must produce a state with norm 0" Δ atol
            success = false
        end
    catch exc
        quiet || @error(
            "$(px)`0.0 * state` must produce a state with norm 0.",
            exception = (exc, catch_abbreviated_backtrace())
        )
        success = false
    end


    if inplace

        has_similar = true
        similar_is_valid = true
        try
            ϕ = similar(state)
            if typeof(ϕ) != typeof(state)
                quiet ||
                    @error "$(px)`similar(state)::$(typeof(ϕ))` must have the same type as `state::$(typeof(state))`"
                success = false
            end
        catch exc
            quiet || @error(
                "$(px)`similar(state)` must be defined.",
                exception = (exc, catch_abbreviated_backtrace())
            )
            has_similar = false
            success = false
        end

        try
            if has_similar
                ϕ = similar(state)
                copyto!(ϕ, state)
                if _check_similar
                    # we only check ϕ after `copyto!`, because just `similar`
                    # might have NaNs in the amplitude, which screws up lots of
                    # the checks
                    if !check_state(
                        ϕ;
                        normalized = false,
                        atol,
                        _check_similar = false,
                        _message_prefix = "On `similar(state)`: ",
                        quiet
                    )
                        quiet || @error("$(px)`similar(state)` must return a valid state")
                        similar_is_valid = false
                        success = false
                    end
                end
                Δ = norm(ϕ - state)
                if Δ > atol
                    quiet ||
                        @error "$(px)`ϕ - state` must have norm 0, where `ϕ = similar(state); copyto!(ϕ, state)`" Δ atol
                    success = false
                end
            end
        catch exc
            quiet || @error(
                "$(px)`copyto!(other, state)` must be defined.",
                exception = (exc, catch_abbreviated_backtrace())
            )
            success = false
        end

        try
            if has_similar && similar_is_valid
                ϕ = similar(state)
                copyto!(ϕ, state)
                state_zero = fill!(ϕ, 0.0)
                if !isa(state_zero, typeof(ϕ))
                    quiet || @error "$(px)`fill!(state, 0.0)` must return the filled state"
                    success = false
                end
                Δ = norm(state_zero)
                if Δ > atol
                    quiet || @error "$(px)`fill!(state, 0.0)` must have norm 0" Δ atol
                    success = false
                end
            end
        catch exc
            quiet || @error(
                "$(px)`fill!(state, c)` must be defined.",
                exception = (exc, catch_abbreviated_backtrace())
            )
            success = false
        end

        try
            if has_similar && similar_is_valid
                ϕ = similar(state)
                copyto!(ϕ, state)
                ϕ = lmul!(1im, ϕ)
                Δ = norm(ϕ - 1im * state)
                if Δ > atol
                    quiet ||
                        @error "$(px)`norm(state)` must have absolute homogeneity: `norm(s * state) = s * norm(state)`" Δ atol
                    success = false
                end
            end
        catch exc
            quiet || @error(
                "$(px)`lmul!(c, state)` for a scalar `c` must be defined.",
                exception = (exc, catch_abbreviated_backtrace())
            )
            success = false
        end

        try
            if has_similar && similar_is_valid
                ϕ = similar(state)
                copyto!(ϕ, state)
                ϕ = lmul!(0.0, ϕ)
                Δ = norm(ϕ)
                if Δ > atol
                    quiet ||
                        @error "$(px)`lmul!(0.0, state)` must produce a state with norm 0" Δ atol
                    success = false
                end
            end
        catch exc
            quiet || @error(
                "$(px)`lmul!(0.0, state)` must produce a state with norm 0.",
                exception = (exc, catch_abbreviated_backtrace())
            )
            success = false
        end

        try
            if has_similar && similar_is_valid
                ϕ = similar(state)
                copyto!(ϕ, state)
                ϕ = axpy!(1im, state, ϕ)
                Δ = norm(ϕ - (state + 1im * state))
                if Δ > atol
                    quiet ||
                        @error "$(px)`axpy!(a, state, ϕ)` must match `ϕ += a * state`" Δ atol
                    success = false
                end
            end
        catch exc
            quiet || @error(
                "$(px)`axpy!(c, state, other)` must be defined.",
                exception = (exc, catch_abbreviated_backtrace())
            )
            success = false
        end

    end  # in-place interface

    if normalized
        try
            η = norm(state)
            Δ = abs(η - 1)
            if Δ > atol
                quiet || @error "$(px)`norm(state)` must be 1, not $η" Δ atol
                success = false
            end
        catch exc
            quiet || @error(
                "$(px)`norm(state)` must be 1.",
                exception = (exc, catch_abbreviated_backtrace())
            )
            success = false
        end
    end

    if supports_vector_interface(state)

        try
            T = eltype(state)
            if !(T isa Type && T <: Number)
                quiet || @error "$(px)`eltype(state)` must return a numeric type, not $T"
                success = false
            end
        catch exc
            quiet || @error(
                "$(px)`eltype(state)` must be defined.",
                exception = (exc, catch_abbreviated_backtrace())
            )
            success = false
        end

        try
            s = size(state)
            if length(s) >= 1 && prod(s) > 0
                val = state[1]
                T = eltype(state)
                if T isa Type && T <: Number && !(val isa T)
                    quiet ||
                        @error "$(px)`state[1]` must return a value of type `eltype(state)=$T`, not $(typeof(val))"
                    success = false
                end
            end
        catch exc
            quiet || @error(
                "$(px)`getindex(state, i)` must be defined.",
                exception = (exc, catch_abbreviated_backtrace())
            )
            success = false
        end

        try
            l = length(state)
            s = size(state)
            if l != prod(s)
                quiet ||
                    @error "$(px)`length(state)` must equal `prod(size(state))`: $l ≠ $(prod(s))"
                success = false
            end
        catch exc
            quiet || @error(
                "$(px)`length(state)` must be defined.",
                exception = (exc, catch_abbreviated_backtrace())
            )
            success = false
        end

        try
            itr = iterate(state)
            s = size(state)
            if isnothing(itr) && prod(s) > 0
                quiet ||
                    @error "$(px)`iterate(state)` must not return `nothing` for a non-empty state"
                success = false
            end
        catch exc
            quiet || @error(
                "$(px)`iterate(state)` must be defined.",
                exception = (exc, catch_abbreviated_backtrace())
            )
            success = false
        end

        try
            st2 = similar(state)
            if !supports_inplace(st2)
                quiet ||
                    @error "$(px)`similar(state)` must return a mutable vector (`supports_inplace` must be `true`), got $(typeof(st2))"
                success = false
            end
            if size(st2) != size(state)
                quiet ||
                    @error "$(px)`similar(state)` must return a vector with the same shape: size $(size(st2)) ≠ $(size(state))"
                success = false
            end
            if eltype(st2) != eltype(state)
                quiet ||
                    @error "$(px)`similar(state)` must return a vector with the same element type: $(eltype(st2)) ≠ $(eltype(state))"
                success = false
            end
        catch exc
            quiet || @error(
                "$(px)`similar(state)` must be defined.",
                exception = (exc, catch_abbreviated_backtrace())
            )
            success = false
        end

        try
            S = (eltype(state) == ComplexF64) ? ComplexF32 : ComplexF64
            st2 = similar(state, S)
            if !supports_inplace(st2)
                quiet ||
                    @error "$(px)`similar(state, $S)` must return a mutable vector (`supports_inplace` must be `true`), got $(typeof(st2))"
                success = false
            end
            if size(st2) != size(state)
                quiet ||
                    @error "$(px)`similar(state, $S)` must return a vector with the same shape: size $(size(st2)) ≠ $(size(state))"
                success = false
            end
            if eltype(st2) != S
                quiet ||
                    @error "$(px)`similar(state, $S)` must return a vector with element type $S, got $(eltype(st2))"
                success = false
            end
        catch exc
            quiet || @error(
                "$(px)`similar(state, ::Type{S})` must be defined.",
                exception = (exc, catch_abbreviated_backtrace())
            )
            success = false
        end

        try
            dims = size(state)
            st2 = similar(state, dims)
            if !supports_inplace(st2)
                quiet ||
                    @error "$(px)`similar(state, dims)` must return a mutable array (`supports_inplace` must be `true`), got $(typeof(st2))"
                success = false
            end
            if size(st2) != dims
                quiet ||
                    @error "$(px)`similar(state, dims)` must return an array with size $dims, got $(size(st2))"
                success = false
            end
            if eltype(st2) != eltype(state)
                quiet ||
                    @error "$(px)`similar(state, dims)` must return an array with the same element type: $(eltype(st2)) ≠ $(eltype(state))"
                success = false
            end
        catch exc
            quiet || @error(
                "$(px)`similar(state, dims::Dims)` must be defined.",
                exception = (exc, catch_abbreviated_backtrace())
            )
            success = false
        end

        try
            S = (eltype(state) == ComplexF64) ? ComplexF32 : ComplexF64
            dims = size(state)
            st2 = similar(state, S, dims)
            if !supports_inplace(st2)
                quiet ||
                    @error "$(px)`similar(state, $S, dims)` must return a mutable array (`supports_inplace` must be `true`), got $(typeof(st2))"
                success = false
            end
            if size(st2) != dims
                quiet ||
                    @error "$(px)`similar(state, $S, dims)` must return an array with size $dims, got $(size(st2))"
                success = false
            end
            if eltype(st2) != S
                quiet ||
                    @error "$(px)`similar(state, $S, dims)` must return an array with element type $S, got $(eltype(st2))"
                success = false
            end
        catch exc
            quiet || @error(
                "$(px)`similar(state, ::Type{S}, dims::Dims)` must be defined.",
                exception = (exc, catch_abbreviated_backtrace())
            )
            success = false
        end

        if supports_inplace(state)

            try
                s = size(state)
                if prod(s) > 0
                    v = state[1]
                    state[1] = v  # write back the same value
                end
            catch exc
                quiet || @error(
                    "$(px)`setindex!(state, v, i)` must be defined for an in-place state.",
                    exception = (exc, catch_abbreviated_backtrace())
                )
                success = false
            end

        end

    end

    return success

end
