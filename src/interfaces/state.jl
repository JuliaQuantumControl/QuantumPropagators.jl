using Test

using LinearAlgebra


"""Check that `state` is a valid element of a Hilbert space.

```julia
@test check_state(state;
                  for_immutable_state=true, for_mutable_state=true,
                  normalized=false, atol=1e-15)
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
"""
function check_state(
    state;
    for_immutable_state=true,
    for_mutable_state=true,
    normalized=false,
    atol=1e-15,
    _check_similar=true  # to avoid infinite recursion
)

    ≈(a, b) = isapprox(a, b; atol)

    success = true

    try
        c = state ⋅ state
        if !(c isa Complex)
            @error "`state ⋅ state` must return a Complex number type, not $(typeof(c))"
            success = false
        end
    catch exc
        @error "the inner product of two states must be a complex number: $exc"
        success = false
    end

    try
        if abs(norm(state) - sqrt(state ⋅ state)) > atol
            @error "`norm(state)` must match √(state⋅state)"
            success = false
        end
    catch exc
        @error "the norm of a state must be defined via the inner product: $exc"
        success = false
    end

    if for_immutable_state

        try
            if norm(state - state) > atol
                @error "`state - state` must have norm 0"
                success = false
            end
            if norm(state + state) > (2 * norm(state) + atol)
                @error "`norm(state + state)` must fulfill the triangle inequality"
                success = false
            end
        catch exc
            @error "`state + state` and `state - state` must be defined: $exc"
            success = false
        end

        try
            ϕ = copy(state)
            if norm(ϕ - state) > atol
                @error "`copy(state) - state` must have norm 0"
                success = false
            end
        catch exc
            @error "copy(state) must be defined: $exc"
            success = false
        end

        try
            ϕ = 0.5 * state
            if abs(norm(ϕ) - 0.5 * norm(state)) > atol
                @error "`norm(state)` must have absolute homogeneity: `norm(s * state) = s * norm(state)`"
                success = false
            end
        catch exc
            @error "`c * state` for a scalar `c` must be defined: $exc"
            success = false
        end

        try
            if norm(0.0 * state) > atol
                @error "`0.0 * state` must produce a state with norm 0"
                success = false
            end
        catch exc
            @error "`0.0 * state` must produce a state with norm 0: $exc"
            success = false
        end

    end

    if for_mutable_state

        if _check_similar
            try
                ϕ = similar(state)
                if !check_state(
                    ϕ;
                    for_immutable_state,
                    for_mutable_state,
                    normalized=false,
                    atol,
                    _check_similar=false
                )
                    @error("`similar(state)` must return a valid state")
                end
            catch exc
                @error "`similar(state)` must be defined: $exc"
                success = false
            end
        end

        try
            ϕ = similar(state)
            copyto!(ϕ, state)
            if for_immutable_state
                if norm(ϕ - state) > atol
                    @error "`copy(state) - state` must have norm 0"
                    success = false
                end
            end
        catch exc
            @error("`copyto!(other, state)` must be defined: $exc")
        end

        try
            ϕ = similar(state)
            copyto!(ϕ, state)
            ϕ = lmul!(1im, ϕ)
            if for_immutable_state
                if norm(ϕ - 1im * state) > atol
                    @error "`norm(state)` must have absolute homogeneity: `norm(s * state) = s * norm(state)`"
                    success = false
                end
            end
        catch exc
            @error("`lmul!(c, state)` for a scalar `c` must be defined: $exc")
        end

        try
            ϕ = similar(state)
            copyto!(ϕ, state)
            ϕ = lmul!(0.0, ϕ)
            if norm(ϕ) > atol
                @error "`lmul!(0.0, state)` must produce a state with norm 0"
                success = false
            end
        catch exc
            @error "`lmul!(0.0, state)` must produce a state with norm 0: $exc"
            success = false
        end

        try
            ϕ = similar(state)
            copyto!(ϕ, state)
            ϕ = axpy!(1im, state, ϕ)
            if for_immutable_state
                if norm(ϕ - (state + 1im * state)) > atol
                    @error "`axpy!(a, state, ϕ)` must match `ϕ += a * state`"
                    success = false
                end
            end
        catch exc
            @error "`axpy!(c, state, other)` must be defined: $exc"
            success = false
        end

    end

    if normalized
        try
            η = norm(state)
            if abs(η - 1) > atol
                @error "`norm(state)` must be 1, not $η"
                success = false
            end
        catch exc
            @error "`norm(state)` must be 1: $exc"
            success = false
        end
    end

    return success

end
