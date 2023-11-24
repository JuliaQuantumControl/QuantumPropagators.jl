using Test

using LinearAlgebra


"""Check that `state` is a valid element of a Hilbert space.

```julia
@test check_state(state;
                  for_immutable_state=true, for_mutable_state=true,
                  normalized=false, atol=1e-15, quiet=false)
```

verifies the following requirements:

* The inner product (`LinearAlgebra.dot`) of two states must return a Complex
  number type.
* The `LinearAlgebra.norm` of `state` must be defined via the inner product.
  This is the *definition* of a Hilbert space, a.k.a a complete inner product
  space or more precisely a Banach space (normed vector space) where the
  norm is induced by an inner product.

If `for_immutable_state`:

* `state + state` and `state - state` must be defined
* `copy(state)` must be defined
* `c * state` for a scalar `c` must be defined
* `norm(state + state)` must fulfill the triangle inequality
* `0.0 * state` must produce a state with norm 0
* `copy(state) - state` must have norm 0
* `norm(state)` must have absolute homogeneity: `norm(s * state) = s *
  norm(state)`

If `for_mutable_state`:

* `similar(state)` must be defined and return a valid state
* `copyto!(other, state)` must be defined
* `LinearAlgebra.lmul!(c, state)` for a scalar `c` must be defined
* `LinearAlgebra.axpy!(c, state, other)` must be defined
* `norm(state)` must fulfill the same general mathematical norm properties as
  with `for_immutable_state`.

If `normalized` (not required by default):

* `LinearAlgebra.norm(state)` must be 1

It is strongly recommended to always support immutable operations (also for
mutable states)

The function returns `true` for a valid state and `false` for an invalid state.
Unless `quiet=true`, it will log an error to indicate which of the conditions
failed.
"""
function check_state(
    state;
    for_immutable_state=true,
    for_mutable_state=true,
    normalized=false,
    atol=1e-14,
    quiet=false,
    _check_similar=true,  # to avoid infinite recursion
    _message_prefix=""  # for recursive calling
)

    ≈(a, b) = isapprox(a, b; atol)
    px = _message_prefix

    success = true

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
            quiet || @error "$(px)`norm(state)=$r1)` must match `√(state⋅state)=$r2` (Δ=$Δ)"
            success = false
        end
    catch exc
        quiet || @error(
            "$(px)the norm of a state must be defined via the inner product.",
            exception = (exc, catch_abbreviated_backtrace())
        )
        success = false
    end

    if for_immutable_state

        try
            Δ = norm(state - state)
            if Δ > atol
                quiet || @error "`$(px)state - state` must have norm 0"
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
            if norm(ϕ - state) > atol
                quiet || @error "$(px)`copy(state) - state` must have norm 0"
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
            if abs(norm(ϕ) - 0.5 * norm(state)) > atol
                quiet ||
                    @error "$(px)`norm(state)` must have absolute homogeneity: `norm(s * state) = s * norm(state)`"
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
            if norm(0.0 * state) > atol
                quiet || @error "$(px)`0.0 * state` must produce a state with norm 0"
                success = false
            end
        catch exc
            quiet || @error(
                "$(px)`0.0 * state` must produce a state with norm 0.",
                exception = (exc, catch_abbreviated_backtrace())
            )
            success = false
        end

    end

    if for_mutable_state

        has_similar = true
        similar_is_valid = true
        try
            ϕ = similar(state)
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
                        for_immutable_state,
                        for_mutable_state,
                        normalized=false,
                        atol,
                        _check_similar=false,
                        _message_prefix="On `similar(state)`: ",
                        quiet
                    )
                        quiet || @error("$(px)`similar(state)` must return a valid state")
                        similar_is_valid = false
                        success = false
                    end
                end
                if for_immutable_state
                    if norm(ϕ - state) > atol
                        quiet ||
                            @error "$(px)`ϕ - state` must have norm 0, where `ϕ = similar(state); copyto!(ϕ, state)`"
                        success = false
                    end
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
                ϕ = lmul!(1im, ϕ)
                if for_immutable_state
                    if norm(ϕ - 1im * state) > atol
                        quiet ||
                            @error "$(px)`norm(state)` must have absolute homogeneity: `norm(s * state) = s * norm(state)`"
                        success = false
                    end
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
                if norm(ϕ) > atol
                    quiet ||
                        @error "$(px)`lmul!(0.0, state)` must produce a state with norm 0"
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
                if for_immutable_state
                    if norm(ϕ - (state + 1im * state)) > atol
                        quiet ||
                            @error "$(px)`axpy!(a, state, ϕ)` must match `ϕ += a * state`"
                        success = false
                    end
                end
            end
        catch exc
            quiet || @error(
                "$(px)`axpy!(c, state, other)` must be defined.",
                exception = (exc, catch_abbreviated_backtrace())
            )
            success = false
        end

    end

    if normalized
        try
            η = norm(state)
            if abs(η - 1) > atol
                quiet || @error "$(px)`norm(state)` must be 1, not $η"
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

    return success

end
