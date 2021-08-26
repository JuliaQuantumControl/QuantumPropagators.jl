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
  `:newton`, and `:expprop`.
"""
function initpropwrk(state, tlist, generator...; method=:auto, kwargs...)
    if method == :auto
        # TODO: determine the best method
        method = :expprop
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


# Default "observable" for storing the propagated state. We want to keep this
# private, as the routine is not safe unless it is the *only*
# observable. If state storage needs to be combined with other observables
# `state->copy(state)` would need to be used.
_store_state(state) = copy(state)
_store_state(state::Vector) = state


"""Propagate a state over an entire time grid.

```julia
propagate(state, genfunc, tlist;
          backwards=false; wrk=nothing, method=:auto, storage=nothing,
          observables=(<store state>, ), hook=nothing)
```

propagates `state` over the time grid in `tlist`, using piecewise-constant
dynamical generators (Hamiltonians or Liouvillians) determined by `genfunc`,
and returns the resulting propagated state. The propagation is performed by
calling [`propstep!`](@ref) for every interval in `tlist`.

For the i'th time interval, `genfunc(tlist, i)` must return the generator for
that time interval. Generally, when approximating a time-continuous dynamical
generator as piecewise-constant on the time grid, it should be evaluated at the
*midpoint* of the interval. A possible exception is the first and last
interval, which may be better evaluated at `tlist[1]` and `tlist[end]` to
ensure exact boundary conditions like control fields that are exactly zero.

In addition to the two positional parameters indicating the time interval,
`genfunc` will also receive the `state` (the input state for the propagation
step), `backwards`, `storage`, `observables`, and `init` as keyword arguments.
These additional parameters may be used for unusual equations of motion beyond
the standard Schrödinger or Liouville-von-Neumann equation, e.g. `state` would
enter the `genfunc` for a Gross–Pitaevskii equation. For standard equations of
motion that do not use the additional parameters, it is best to capture the
keyword arguments to `genfunc` with a definition like

```julia
genfunc(tlist, i; kwargs...) = ...
```

The propagation method is determined by `wrk`, see [`initpropwrk`](@ref) and
[`propstep!](@ref): If `wrk` is not given, it will be created internally for
the given `method` with [`initpropwrk`](@ref); for the default `method=:auto`,
[`initpropwrk`](@ref) will attempt to choose the most appropriate method for
the given parameters, using the generator returned by `genfunc` with `i=1` and
the keyword argument `init=true`.

In general, there is no requirement that `tlist` has a constant time step,
although some propagation methods (most notably [`cheby!`](@ref)) only support
a uniform time grid.

If `storage` is given as an Array, it will be filled with data determined by
the `observables`. The default "observable" results in the propagated states at
every point in time being stored.
The `storage` array should be created with [`init_storage`](@ref). See its
documentation for details.

The `storage` parameter may also be given as `true`, and a new storage array
will be created internally with [`init_storage`](@ref) and returned instead of
the propagated state.

If `backwards` is `true`, the input state is assumed to be at time
`tlist[end]`, and the propagation progresses backwards in time (with a negative
time step `dt`). If `storage` is given, it will be filled back-to-front during
the backwards propagation.

If `hook` is given as a callable, it will be called after each propagation
step, as `hook(state, generator, tlist, i, wrk, observables)` where `i` is the
index of the time interval on `tlist` covered by the propagation step (0 for
the initial state, respectives `lastindex(tlist)` for the backward
propagation).  The `hook` is called before calculating any observables. Example
usage includes writing data to file, or modyfing `state`, e.g. removing
amplitude from the lowest and highest level to mitigate "truncation error".

The `propagate` routine returns the propagated state at `tlist[end]`,
respectively `tlist[1]` if `backwards=true`, or a storage array with the
stored states / observable data if `storage=true`.
"""
function propagate(state, genfunc, tlist;
                   backwards=false,
                   wrk=nothing,
                   method=:auto,
                   storage=nothing,
                   observables=(_store_state, ),
                   hook=nothing,
                  )
    return_storage = false
    if storage === true
        storage = init_storage(state, tlist, observables)
        return_storage = true
    end
    state = copy(state)
    if wrk == nothing
        G = genfunc(tlist, 1; state=state, backwards=backwards,
                    storage=storage, observables=observables, init=true)
        wrk = initpropwrk(state, tlist, G; method=method)
    end
    intervals = enumerate(tlist[2:end])
    if backwards
        intervals = Iterators.reverse(intervals)
        if hook ≠ nothing
            hook(state, generator, tlist, lastindex(tlist), wrk, observables)
        end
        if storage ≠ nothing
            write_to_storage!(storage, lastindex(tlist), state, observables)
        end
    else
        if hook ≠ nothing
            hook(state, generator, tlist, 0, wrk, observables)
        end
        if storage ≠ nothing
            write_to_storage!(storage, 1, state, observables)
        end
    end
    # TODO: optional progress meter
    for (i, t_end) in intervals
        dt = t_end - tlist[i]
        generator = genfunc(tlist, i; state=state, backwards=backwards,
                            storage=storage, observables=observables,
                            init=false)
        propstep!(state, generator, (backwards ? -dt : dt), wrk)
        if storage ≠ nothing
            write_to_storage!(storage, i + (backwards ? 0 : 1), state,
                              observables)
        end
        if hook ≠ nothing
            hook(state, generator, tlist, i, wrk, observables)
        end
    end
    if return_storage
        return storage
    else
        return state
    end
end
