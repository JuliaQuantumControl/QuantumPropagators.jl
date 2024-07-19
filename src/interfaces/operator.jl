using Test

using LinearAlgebra
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

* `op` must not be time-dependent: [`evaluate(op, tlist, 1) ≡ op`](@ref evaluate)
* `op` must not contain any controls: [`length(get_controls(op)) == 0`](@ref get_controls)
* `op * state` must be defined
* The [`QuantumPropagators.Interfaces.supports_inplace`](@ref) method must be
  defined for `op`. If it returns `true`, it must be possible to evaluate a
  generator in-place into the existing `op`. See [`check_generator`](@ref).

If [`QuantumPropagators.Interfaces.supports_inplace(state)`](@ref):

* The 3-argument `LinearAlgebra.mul!` must apply `op` to the given `state`
* The 5-argument `LinearAlgebra.mul!` must apply `op` to the given `state`
* `LinearAlgebra.mul!` must match `*`, if applicable
* `LinearAlgebra.mul!` must return the resulting state

If `for_expval` (typically required for optimal control):

* `LinearAlgebra.dot(state, op, state)` must return return a number
* `dot(state, op, state)` must match `dot(state, op * state)`, if applicable

The function returns `true` for a valid operator and `false` for an invalid
operator. Unless `quiet=true`, it will log an error to indicate which of the
conditions failed.
"""
function check_operator(
    op;
    state,
    tlist=[0.0, 1.0],
    for_expval=true,
    atol=1e-14,
    quiet=false,
    _message_prefix=""  # for recursive calling
)

    ≈(a, b) = isapprox(a, b; atol)

    px = _message_prefix
    success = true

    Ψ = state

    @assert check_state(state; atol, quiet=true)
    @assert tlist isa Vector{Float64}

    try
        supports_inplace(op)
    catch exc
        quiet || @error(
            "$(px)The `QuantumPropagators.Interfaces.supports_inplace` method must be defined for `op`.",
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
            if norm(ϕ - op * Ψ) > atol
                quiet || @error "$(px)`mul!(ϕ, op, state)` must match `op * state`"
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
            if norm(ϕ - (0.5 * ϕ0 + 0.5 * (op * Ψ))) > atol
                quiet ||
                    @error "$(px)`mul!(ϕ, op, state, α, β)` must match β*ϕ + α*op*state"
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
            if abs(dot(Ψ, op, Ψ) - dot(Ψ, op * Ψ)) > atol
                quiet ||
                    @error "$(px)`dot(state, op, state)` must match `dot(state, op * state)`"
                success = false
            end
            if supports_inplace(state)
                ϕ = similar(Ψ)
                mul!(ϕ, op, Ψ)
                if abs(dot(Ψ, op, Ψ) - dot(Ψ, ϕ)) > atol
                    quiet ||
                        @error "$(px)`dot(state, op, state)` must match `dot(state, op * state)`"
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

    return success

end
