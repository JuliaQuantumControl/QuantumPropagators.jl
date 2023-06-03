import ..AbstractPropagator
import ..PiecewisePropagator
import ..init_prop
import ..reinit_prop!
import ..prop_step!
import ..set_state!
import ..set_t!

"""Check that the given `propagator` implements the required interface.

```julia
@test check_propagator(propagator; atol=1e-14)
```

verifies that the `propagator` matches the interface described for an
[`AbstractPropagator`](@ref). The `propagator` must have been freshly
initialized with [`init_prop`](@ref).

* `propagator` must have the properties `state`, `tlist`, `t`, `backward`, and
  `inplace`
* `propagator.state` must be a valid state (see [`check_state`](@ref)), with
  support for in-place operations (`for_mutable_state=true`) if
  `propagator.inplace` is true.
* `propagator.t` must be the first or last element of `propagator.tlist`,
  depending on `propagator.backward`
* [`prop_step!(propagator)`](@ref prop_step!) must be defined and return a
  valid state until the time grid is exhausted
* For an in-place propagator, the state returned by [`prop_step!`](@ref) must
  be the `propagator.state` object
* For a not-in-place propagator, the state returned by [`prop_step!`](@ref)
  must be a new object
* [`prop_step!`](@ref) must advance `propagator.t` forward or backward one step
  on the time grid
* [`prop_step!`](@ref) must return `nothing` when going beyond the time grid
* [`set_t!(propagator, t)`](@ref set_t!) must be defined and set `propagator.t`
* [`set_state!(propagator, state)`](@ref set_state!) must be defined and set
  `propagator.state`.
* [`set_state!(propagator, state)`](@ref set_state!) for an in-place propagator
    must overwrite `propagator.state` in-place.
* In a [`PiecewisePropagator`](@ref), `propagator.parameters` must be a dict
  mapping controls to a vector of values, one for each interval on
  `propagator.tlist`
* [`reinit_prop!`](@ref) must be defined and re-initialize the propagator
* [`reinit_prop!(propagator, state)`](@ref reinit_prop!) must be idempotent.
  That is, repeated calls to [`reinit_prop!`](@ref) leave the `propagator`
  unchanged.
"""
function check_propagator(propagator; atol=1e-14)

    local state, tlist, t, parameters, backward, inplace

    success = true

    try
        state = propagator.state
        tlist = propagator.tlist
        t = propagator.t
        parameters = propagator.parameters
        backward = propagator.backward
        inplace = propagator.inplace
    catch exc
        @error "`propagator` does not have the required properties: $exc"
        success = false
    end

    success || return false  # no point in going on

    try
        if !(check_state(state; for_mutable_state=inplace, for_immutable_state=true, atol))
            @error "`propagator.state` is not a valid state"
            success = false
        end
    catch exc
        @error "Cannot verify `propagator.state`: $exc"
        success = false
    end

    try
        if backward
            if t ≠ tlist[end]
                @error "propagator.t ≠ propagator.tlist[end]"
                success = false
            end
        else
            if t ≠ tlist[begin]
                @error "propagator.t ≠ propagator.tlist[begin]"
                success = false
            end
        end
    catch exc
        @error "Cannot verify `propagator.t`: $exc"
        success = false
    end

    Ψ₀_ref = state
    Ψ₀ = copy(state)
    Ψ₁ = similar(Ψ₀)

    try
        Ψ₁ = prop_step!(propagator)
        if !(check_state(Ψ₁; for_mutable_state=false, for_immutable_state=false, atol))
            @error "prop_step! must return a valid state until time grid is exhausted"
            success = false
        end
        if inplace
            if Ψ₁ ≢ state
                @error "For an in-place propagator, the state returned by `prop_step!` must be the `propagator.state` object"
                success = false
            end
        else
            if Ψ₁ ≡ state
                @error "For a not-in-place propagator, the state returned by `prop_step!` must be a new object"
                success = false
            end
        end
        if backward
            Δ = propagator.t - tlist[end-1]
        else
            Δ = propagator.t - tlist[begin+1]
        end
        if abs(Δ) > atol
            @error "`prop_step!` must advance `propagator.t` forward or backward one step on the time grid (Δ=$Δ)"
            success = false
        end
    catch exc
        @error "`prop_step!(propagator)` must be defined : $exc"
        success = false
    end

    try
        t = tlist[end]
        if backward
            t = tlist[begin]
        end
        set_t!(propagator, t)
        if propagator.t ≠ t
            @error "`set_t!(propagator, t)` must set propagator.t=$t, not $(propagator.t)"
            success = false
        end
    catch exc
        @error "`set_t!(propagator)` must be defined : $exc"
        success = false
    end

    try
        t = propagator.t
        res = prop_step!(propagator)
        if !isnothing(res)
            success = false
            if backward
                @error "`prop_step!(propagator)` at initial t=$t must return `nothing`"
            else
                @error "`prop_step!(propagator)` at final t=$t must return `nothing`"
            end
        end
    catch exc
        @error "Failed to run `prop_step!(propagator)` at t=$(propagator.t): $exc"
        success = false
    end

    try
        set_state!(propagator, Ψ₀)
        if norm(propagator.state - Ψ₀) > atol
            @error "`set_state!(propagator, state)` must set `propagator.state`"
        end
        if inplace
            if propagator.state ≢ Ψ₀_ref
                @error "`set_state!(propagator, state)` for an in-place propagator must overwrite `propagator.state` in-place."
            end
        end
    catch exc
        @error "`set_state!(propagator, state)` must be defined : $exc"
        success = false
    end

    if propagator isa PiecewisePropagator
        try
            for (control, ampl) ∈ propagator.parameters
                if length(ampl) ≠ length(propagator.tlist) - 1
                    @error "In a PiecewisePropagator, `propagator.parameters` must be a dict mapping controls to a vector of values, one for each interval on `propagator.tlist`"
                    success = false
                    break
                end
            end
        catch exc
            @error "In a PiecewisePropagator, `propagator.parameters` must be a dict mapping controls to a vector of values, one for each interval on `propagator.tlist`: $exc"
            success = false
        end
    end

    try
        reinit_prop!(propagator, Ψ₀)
        t = backward ? propagator.tlist[end] : propagator.tlist[begin]
        if abs(propagator.t - t) > atol
            @error "`reinit_prop!(propagator, state)` must reset `propagator.t`"
            success = false
        end
        if norm(propagator.state - Ψ₀) > atol
            @error "`reinit_prop!(propagator, state)` must reset `propagator.state`"
            success = false
        end
        if propagator.backward ≠ backward
            @error "`reinit_prop!(propagator, state)` must keep `propagator.backward`"
            success = false
        end
        if propagator.inplace ≠ inplace
            @error "`reinit_prop!(propagator, state)` must keep `propagator.inplace`"
            success = false
        end
        prop_step!(propagator)
        if norm(propagator.state - Ψ₁) > atol  # Ψ₁ from very first prop_step!
            # Without keyword arguments to reinit_prop!, the propagator should
            # be in the exact original state and thus reproduce the exact same
            # propagation
            @error "`reinit_prop!(propagator, state)` must be idempotent"
            success = false
        end
        if inplace
            reinit_prop!(propagator, Ψ₀)
        else
            @assert norm(Ψ₀_ref - Ψ₀) < atol
            reinit_prop!(propagator, Ψ₀_ref)
            if propagator.state ≢ Ψ₀_ref
                @error "`reinit_prop!(propagator, state)` for in-place propagator must be in-place"
                success = false
            end
        end
    catch exc
        @error "`reinit_prop!` must be defined and re-initialize the propagator: $exc"
        success = false
    end

    return success

end
