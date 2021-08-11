module QuantumPropagators

include("./cheby.jl")
include("./newton.jl")
include("./expprop.jl")

export cheby_coeffs, cheby_coeffs!, ChebyWrk, cheby!
export NewtonWrk, newton!
export ExpPropWrk, expprop!
export initpropwrk, propstep!, propagate


"""Initialize a workspace for propagation.

```julia
    wrk = initpropwrk(state, tlist, generator...; method=:auto, kwargs...)
```

The resulting `wrk` can be passed to [`propagate`](@ref) or
[`propstep!`](@ref).

# Arguments

* `state`: An exemplary state for the propagation (e.g., the initial state)
* `tlist`: The time grid over which [`propagate`](@ref) will be called. Must
  include at least to points in order to determine the propagation time step
  to prepare. If the propagation will be over a `tlist` with a variable `dt`,
  the full `tlist` must be passed here.
* `generator`: An exemplary (non-time-dependent) dynamical generator. For full
  generality (if `method=:cheby`), the given `generator` should have a spectral
  range sufficiently large to encompass the entire propagation. If given
  multiple times, a spectral envelope enclosing all the generators will be
  determined automatically. In this case, you should pass the generators with
  the extremal values of all the controls.
* `method`: The propagation method to use. The default value of `:auto`
  attempts to choose the best method available, based on the properties of the
  given `state`, `tlist`, and `generator`. Alternative values are `:cheby` and
  `:newton`.
"""
function initpropwrk(state, tlist, generator...; method=:auto, kwargs...)
    if method == :auto
        # TODO: determine the best method
        error("not implemented")
    end
    if method == :cheby
        # TODO Find the spectral envelope of all generators
        # TODO ensure that `tlist` has a constant dt
        Δ = 0
        E_min = 0
        dt = tlist[2] - tlist[1]
        return ChebyWrk(state, Δ, E_min, dt; kwargs...)
    elseif method == :newton
        return NewtonWrk(state; kwargs...)
    elseif method == :expprop
        return ExpPropWrk(state; kwargs...)
    else
        throw(ArgumentError("Unknown method $method"))
    end
end


"""Perform a single propagation step in-place.

```julia
    propstep!(state, generator, dt, wrk;, kwargs...)
```

The propagation method is determined by `wrk`, see [`initpropwrk`](@ref).
The `kwargs` are forwarded to the underlying method
"""
function propstep!(state, generator, dt, wrk::ChebyWrk; kwargs...)
    return cheby!(state, generator, dt, wrk; kwargs...)
end

function propstep!(state, generator, dt, wrk::NewtonWrk; kwargs...)
    return newton!(state, generator, dt, wrk; kwargs...)
end

function propstep!(state, generator, dt, wrk::ExpPropWrk; kwargs...)
    return expprop!(state, generator, dt, wrk; kwargs...)
end


"""Propagate a state over an entire time grid.

```julia
    propagate(state, genfunc, tlist; wrk=nothing, storage=nothing)
```

Propagates `state` over the `tlist` time grid, using `genfunc(t)` as the
dynamical generator with `t` being the midpoint of each time interval. Return
the propagated state at `tlist[end]`.

The propagation method is determined by `wrk`, see [`initpropwrk`](@ref). If
`wrk` is not given, `propagate` will attempt to choose the most appropriate
method for the given parameters.

If `storage` is given, the propagated states at each point in `tlist` will be
written to it.
"""
function propagate(state, genfunc, tlist; wrk=nothing, storage=nothing)
    state = copy(state)
    if wrk == nothing
        wrk = initipropwrk(state, tlist, genfunc(0))
    end
    # TODO: optional progress meter
    for (interval, t_end) in enumerate(tlist[2:end])
        dt = t_end - tlist[interval-1]
        t = t_end - 0.5*dt
        generator = genfunc(t) # TODO might be better to tell it `interval`
        propstep!(state, generator, dt, wrk)
        # TODO: copy state into storage if requested
    end
    return state
end

end
