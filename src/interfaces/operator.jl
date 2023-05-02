using Test

using LinearAlgebra
using ..Controls: get_controls, evaluate


"""Check that `op` is a valid operator that can be applied to `state`.

```julia
@test check_operator(op; state, tlist=[0.0, 1.0],
                     for_mutable_state=true, for_immutable_state=true,
                     for_expval=true, atol=1e-15)
```

verifies the given `op` relative to `state`. The `state` must pass
[`check_state`](@ref).

An "operator" is any object that [`evaluate`](@ref) returns when evaluating a
time-dependent dynamic generator. The specific requirements for `op` are:

* `op` must not be time-dependent: [`evaluate(op, tlist, 1) ≡ op`](@ref evaluate)
* `op` must not contain any controls: [`length(get_controls(op)) == 0`](@ref get_controls)

If `for_immutable_state` (e.g., for use in propagators with `inplace=false`):

* `op * state` must be defined

If `for_mutable_state` (e.g., for use in propagators with `inplace=true`):

* The 3-argument `LinearAlgebra.mul!` must apply `op` to the given `state`
* The 5-argument `LinearAlgebra.mul!` must apply `op` to the given `state`
* `LinearAlgebra.mul!` must match `*`, if applicable
* `LinearAlgebra.mul!` must return the resulting state

If `for_expval` (typically required for optimal control):

* `LinearAlgebra.dot(state, op, state)` must return return a number
* `dot(state, op, state)` must match `dot(state, op * state)`, if applicable
"""
function check_operator(
    op;
    state,
    tlist=[0.0, 1.0],
    for_mutable_state=true,
    for_immutable_state=true,
    for_expval=true,
    atol=1e-15
)

    ≈(a, b) = isapprox(a, b; atol)

    success = true

    Ψ = state

    @assert check_state(state; for_mutable_state, for_immutable_state, atol)
    @assert tlist isa Vector{Float64}

    try
        H = evaluate(op, tlist, 1)
        if H ≢ op
            @error "`evaluate(op, tlist, 1) must return operator ≡ op"
            success = false
        end
    catch exc
        @error "op must not be time-dependent: $exc"
        success = false
    end

    try
        controls = get_controls(op)
        if length(controls) ≠ 0
            @error "get_controls(op) must return an empty tuple, not $controls"
            success = false
        end
    catch exc
        @error "op must not contain any controls: $exc"
        success = false
    end

    if for_immutable_state

        try
            ϕ = op * Ψ
            if !(ϕ isa typeof(Ψ))
                @error "`op * state` must return an object of the same type as `state`, not $typeof(ϕ)"
                success = false
            end
        catch exc
            @error "`op * state` must be defined: $exc"
            success = false
        end

    end

    if for_mutable_state

        try
            ϕ = similar(Ψ)
            if mul!(ϕ, op, Ψ) ≢ ϕ
                @error "`mul!(ϕ, op, state)` must return the resulting ϕ"
                success = false
            end
            if for_immutable_state
                if norm(ϕ - op * Ψ) > atol
                    @error "`mul!(ϕ, op, state)` must match `op * state`"
                    success = false
                end
            end
        catch exc
            @error "The 3-argument `mul!` must apply `op` to the given `state`: $exc"
            success = false
        end

        try
            ϕ = similar(Ψ)
            ϕ0 = similar(Ψ)
            copyto!(ϕ, Ψ)
            copyto!(ϕ0, Ψ)
            if mul!(ϕ, op, Ψ, 0.5, 0.5) ≢ ϕ
                @error "`mul!(ϕ, op, state, α, β)` must return the resulting ϕ"
                success = false
            end
            if for_immutable_state
                if norm(ϕ - (0.5 * ϕ0 + 0.5 * (op * Ψ))) > atol
                    @error "`mul!(ϕ, op, state, α, β)` must match β*ϕ + α*op*state"
                    success = false
                end
            end
        catch exc
            @error "The 5-argument `mul!` must apply `op` to the given `state`: $exc"
            success = false
        end

    end

    if for_expval

        try
            val = dot(Ψ, op, Ψ)
            if !(val isa Number)
                @error "`dot(state, op, state)` must return a number, not $typeof(val)"
                success = false
            end
            if for_immutable_state
                if abs(dot(Ψ, op, Ψ) - dot(Ψ, op * Ψ)) > atol
                    @error "`dot(state, op, state)` must match `dot(state, op * state)`"
                    success = false
                end
            end
            if for_mutable_state
                ϕ = similar(Ψ)
                mul!(ϕ, op, Ψ)
                if abs(dot(Ψ, op, Ψ) - dot(Ψ, ϕ)) > atol
                    @error "`dot(state, op, state)` must match `dot(state, op * state)`"
                    success = false
                end
            end
        catch exc
            @error "`dot(state, op, state)` must return return a number: $exc"
            success = false
        end

    end

    return success

end
