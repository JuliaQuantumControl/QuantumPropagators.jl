import ..AbstractPropagator
import ..PiecewisePropagator
import ..init_prop
import ..reinit_prop!
import ..prop_step!
import ..set_state!
import ..set_t!


"""Check that the given `propagator` implements the required interface.

```julia
@test check_propagator(propagator; atol=1e-14, quiet=false)
```

verifies that the `propagator` matches the interface described for an
[`AbstractPropagator`](@ref). The `propagator` must have been freshly
initialized with [`init_prop`](@ref).

* `propagator` must have the properties `state`, `tlist`, `t`, `parameters`,
  `backward`, and `inplace`
* `propagator.state` must be a valid state (see [`check_state`](@ref))
* If `propagator.inplace` is true, [`supports_inplace`](@ref) for
  `propagator.state` must also be true
* `propagator.tlist` must be monotonically increasing.
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
* [`set_state!`](@ref) must return the set `propagator.state`
* In a [`PiecewisePropagator`](@ref), `propagator.parameters` must be a dict
  mapping controls to a vector of values, one for each interval on
  `propagator.tlist`
* [`reinit_prop!`](@ref) must be defined and re-initialize the propagator
* [`reinit_prop!(propagator, state)`](@ref reinit_prop!) must be idempotent.
  That is, repeated calls to [`reinit_prop!`](@ref) leave the `propagator`
  unchanged.

The function returns `true` for a valid propagator and `false` for an invalid
propagator. Unless `quiet=true`, it will log an error to indicate which of the
conditions failed.
"""
function check_propagator(
    propagator;
    atol=1e-14,
    quiet=false,
    _message_prefix=""  # for recursive calling
)

    local state, tlist, t, parameters, backward, inplace

    px = _message_prefix
    success = true

    try
        state = propagator.state
        tlist = propagator.tlist
        t = propagator.t
        parameters = propagator.parameters
        backward = propagator.backward
        inplace = propagator.inplace
    catch exc
        quiet || @error(
            "$(px)`propagator` does not have the required properties.",
            exception = (exc, catch_abbreviated_backtrace())
        )
        success = false
    end

    success || return false  # no point in going on

    valid_tlist =
        check_tlist(propagator.tlist; quiet, _message_prefix="On `propagator.tlist`: ")
    if !valid_tlist
        quiet || @error "$(px)`propagator.tlist` is not a valid time grid"
        success = false
    end

    try
        valid_state = check_state(state; atol, _message_prefix="On `propagator.state`: ")
        if !valid_state
            quiet || @error "$(px)`propagator.state` is not a valid state"
            success = false
        end
        if inplace
            if !supports_inplace(state)
                quiet ||
                    @error "$(px)`supports_inplace(propagator.state)` is false even for an in-place `propagator`"
                success = false
            end
        end
    catch exc
        quiet || @error(
            "$(px)Cannot verify `propagator.state`.",
            exception = (exc, catch_abbreviated_backtrace())
        )
        success = false
    end

    try
        _t = tlist[begin]
        for t in tlist[begin+1:end]
            if t <= _t
                quiet || @error "$(px)`propagator.tlist` must be monotonically increasing"
                success = false
            end
            _t = t
        end
    catch exc
        quiet || @error(
            "$(px)Cannot verify `propagator.tlist`.",
            exception = (exc, catch_abbreviated_backtrace())
        )
        success = false
    end

    try
        if backward
            if t ≠ tlist[end]
                quiet || @error "$(px)propagator.t ≠ propagator.tlist[end]"
                success = false
            end
        else
            if t ≠ tlist[begin]
                quiet || @error "$(px)propagator.t ≠ propagator.tlist[begin]"
                success = false
            end
        end
    catch exc
        quiet || @error(
            "$(px)Cannot verify `propagator.t`.",
            exception = (exc, catch_abbreviated_backtrace())
        )
        success = false
    end

    Ψ₀_ref = state
    Ψ₀ = copy(state)
    Ψ₁ = similar(Ψ₀)

    try
        Ψ₁ = prop_step!(propagator)
        valid_state =
            check_state(Ψ₁; atol, _message_prefix="On `Ψ₁=prop_step!(propagator)`: ")
        if !valid_state
            quiet ||
                @error "$(px)prop_step! must return a valid state until time grid is exhausted"
            success = false
        end
        if inplace
            if Ψ₁ ≢ state
                quiet ||
                    @error "$(px)For an in-place propagator, the state returned by `prop_step!` must be the `propagator.state` object"
                success = false
            end
        else
            if Ψ₁ ≡ state
                quiet ||
                    @error "$(px)For a not-in-place propagator, the state returned by `prop_step!` must be a new object"
                success = false
            end
        end
        if backward
            Δ = propagator.t - tlist[end-1]
        else
            Δ = propagator.t - tlist[begin+1]
        end
        if abs(Δ) > atol
            quiet ||
                @error "$(px)`prop_step!` must advance `propagator.t` forward or backward one step on the time grid" Δ atol
            success = false
        end
    catch exc
        quiet || @error(
            "$(px)`prop_step!(propagator)` must be defined.",
            exception = (exc, catch_abbreviated_backtrace())
        )
        success = false
    end

    try
        t = tlist[end]
        if backward
            t = tlist[begin]
        end
        set_t!(propagator, t)
        if propagator.t ≠ t
            quiet ||
                @error "$(px)`set_t!(propagator, t)` must set propagator.t=$t, not $(propagator.t)"
            success = false
        end
    catch exc
        quiet || @error(
            "$(px)`set_t!(propagator)` must be defined.",
            exception = (exc, catch_abbreviated_backtrace())
        )
        success = false
    end

    try
        t = propagator.t
        res = prop_step!(propagator)
        if !isnothing(res)
            success = false
            if backward
                quiet ||
                    @error "$(px)`prop_step!(propagator)` at initial t=$t must return `nothing`"
                success = false
            else
                quiet ||
                    @error "$(px)`prop_step!(propagator)` at final t=$t must return `nothing`"
                success = false
            end
        end
    catch exc
        quiet || @error(
            "$(px)Failed to run `prop_step!(propagator)` at t=$(propagator.t).",
            exception = (exc, catch_abbreviated_backtrace())
        )
        success = false
    end

    try
        Ψ = set_state!(propagator, Ψ₀)
        Δ = norm(propagator.state - Ψ₀)
        if Δ > atol
            quiet ||
                @error "$(px)`set_state!(propagator, state)` must set `propagator.state`" Δ atol
            success = false
        end
        if inplace
            if propagator.state ≢ Ψ₀_ref
                quiet ||
                    @error "$(px)`set_state!(propagator, state)` for an in-place propagator must overwrite `propagator.state` in-place."
                success = false
            end
        end
        if propagator.state ≢ Ψ
            quiet ||
                @error "$(px)`set_state!(propagator, state)` must return `propagator.state`."
            success = false
        end
    catch exc
        quiet || @error(
            "$(px)`set_state!(propagator, state)` must be defined.",
            exception = (exc, catch_abbreviated_backtrace())
        )
        success = false
    end

    if propagator isa PiecewisePropagator
        try
            for (control, ampl) ∈ propagator.parameters
                if length(ampl) ≠ length(propagator.tlist) - 1
                    quiet ||
                        @error "$(px)In a PiecewisePropagator, `propagator.parameters` must be a dict mapping controls to a vector of values, one for each interval on `propagator.tlist`"
                    success = false
                    break
                end
            end
        catch exc
            quiet || @error(
                "$(px)In a PiecewisePropagator, `propagator.parameters` must be a dict mapping controls to a vector of values, one for each interval on `propagator.tlist`.",
                exception = (exc, catch_abbreviated_backtrace())
            )
            success = false
        end
    end

    try
        reinit_prop!(propagator, Ψ₀)
        t = backward ? propagator.tlist[end] : propagator.tlist[begin]
        Δ = abs(propagator.t - t)
        if Δ > atol
            quiet ||
                @error "$(px)`reinit_prop!(propagator, state)` must reset `propagator.t`" Δ atol
            success = false
        end
        Δ = norm(propagator.state - Ψ₀)
        if Δ > atol
            quiet ||
                @error "$(px)`reinit_prop!(propagator, state)` must reset `propagator.state`" Δ atol
            success = false
        end
        if propagator.backward ≠ backward
            quiet ||
                @error "$(px)`reinit_prop!(propagator, state)` must keep `propagator.backward`"
            success = false
        end
        if propagator.inplace ≠ inplace
            quiet ||
                @error "$(px)`reinit_prop!(propagator, state)` must keep `propagator.inplace`"
            success = false
        end
        prop_step!(propagator)
        Δ = norm(propagator.state - Ψ₁)
        if Δ > atol  # Ψ₁ from very first prop_step!
            # Without keyword arguments to reinit_prop!, the propagator should
            # be in the exact original state and thus reproduce the exact same
            # propagation
            quiet ||
                @error "$(px)`reinit_prop!(propagator, state)` must be idempotent" Δ atol
            success = false
        end
        if inplace
            reinit_prop!(propagator, Ψ₀)
        else
            @assert norm(Ψ₀_ref - Ψ₀) < atol
            reinit_prop!(propagator, Ψ₀_ref)
            if propagator.state ≢ Ψ₀_ref
                quiet ||
                    @error "$(px)`reinit_prop!(propagator, state)` for in-place propagator must be in-place"
                success = false
            end
        end
    catch exc
        quiet || @error(
            "$(px)`reinit_prop!` must be defined and re-initialize the propagator.",
            exception = (exc, catch_abbreviated_backtrace())
        )
        success = false
    end

    return success

end
