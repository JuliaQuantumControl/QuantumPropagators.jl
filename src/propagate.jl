using LinearAlgebra

"""Initialize a workspace for propagation.

```julia
wrk = initpropwrk(state, tlist, method=:auto, generator...; kwargs...)
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
* `specrad_method`: for `method=:cheby`, method to use for estimating the
   spectral radius, see [`specrange`](@ref). Defaults to `:auto`.
* `tolerance`: for `method=:cheby`, a tolerance factor for the estimated
  spectral radius. That is, Chebychev coefficients will be calculated for a
  spectral radius increased by the `tolerance` factor compared to the specral
  radius estimated for the `generator`.

All other `kwargs` are filtered and passed to the contructor for returned
workspace, e.g. `limit` for `method=:cheby` or `m_max` for `method=:newton`.
For `method=:cheby`, they additionally passed to [`specrange`](@ref).
"""
function initpropwrk(state, tlist, method::Val{:auto}, generator...; kwargs...)
    is_small = false
    if length(generator) == 0
        try
            is_small = length(state) <= 8
        catch exception
            if isa(exception, MethodError)
                @warn "Cannot determine problem dimension"
            else
                rethrow()
            end
        end
        method = is_small ? Val(:expprop) : Val(:newton)
    else # we have at least one generator
        try
            is_small = length(state) <= 8
        catch
            try
                is_small = all(size(G)[2] <= 8 for G in generator)
            catch exception
                if isa(exception, MethodError)
                    @warn "Cannot determine problem dimension"
                else
                    rethrow()
                end
            end
        end
        if is_small
            method = Val(:expprop)
        else
            method = Val(:newton)
            try
                if all(ishermitian(G) for G in generator)
                    method = Val(:cheby)
                end
            catch exception
                # if we can't determine G is Hermitian, we don't do anything
                if isa(exception, MethodError)
                    @warn("Cannot determine if generator is Hermitian.")
                else
                    rethrow()
                end
            end
        end
    end
    return initpropwrk(state, tlist, method, generator...; kwargs)
end

function initpropwrk(state, tlist, method, generator...; kwargs...)
    throw(TypeError("`method` must be a known symbol, not $method"))
end


function initpropwrk(state, tlist; kwargs...)
    return initpropwrk(state, tlist, Val(:auto); kwargs...)
end


function initpropwrk(state, tlist, method::Symbol, generator...; kwargs...)
    return initpropwrk(state, tlist, Val(method), generator...; kwargs...)
end


function initpropwrk(state, tlist, method::Val{:cheby}, generator...;
                     specrad_method=Val(:auto), tolerance=0.01, kwargs...)
    if length(generator) == 0
        throw(ArgumentError(
            "The :cheby method requires at least on `generator`"
        ))
    end
    # find a spectral envelope for all generators
    E_min, E_max = specrange(generator[1], specrad_method; kwargs...)
    for G in generator[2:end]
        _E_min, _E_max = specrange(G, specrad_method; kwargs...)
        E_min = (_E_min < E_min) ? _E_min : E_min
        E_max = (_E_max > E_max) ? _E_max : E_max
    end
    Δ = E_max - E_min
    δ = 0.01 * Δ
    E_min = E_min - δ/2
    Δ = Δ + δ
    dt = tlist[2] - tlist[1]
    allowed_kwargs = Set((:limit, ))
    filtered_kwargs = filter(p->p.first in allowed_kwargs, kwargs)
    return ChebyWrk(state, Δ, E_min, dt; filtered_kwargs...)
end


function initpropwrk(state, tlist, method::Val{:newton}, generator...; kwargs...)
    allowed_kwargs = Set((:m_max, ))
    filtered_kwargs = filter(p->p.first in allowed_kwargs, kwargs)
    return NewtonWrk(state; filtered_kwargs...)
end


function initpropwrk(state, tlist, method::Val{:expprop}, generator...; kwargs...)
    # ExpPropWrk has no kwargs
    return ExpPropWrk(state)
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
state_out = propagate(
    state, genfunc, tlist, method=:auto;
    backwards=false; storage=nothing, observables=(<store state>, ),
    hook=nothing, kwargs...)
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

For valid propagation `method`s, see [`initpropwrk`](@ref).

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
the propagated state:

```julia
data = propagate(
    state, genfunc, tlist; method=:auto
    backwards=false; storage=true, observables=observables,
    hook=nothing, kwargs...)
```

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
function propagate(state, genfunc, tlist; method=Val(:auto), kwargs...)
    return propagate(state, genfunc, tlist, method; kwargs...)
end

function propagate(state, genfunc, tlist, method::Val; kwargs...)
    # TODO: document what happens with kwargs
    backwards = get(kwargs, :backwards, false)
    storage = get(kwargs, :storage, nothing)
    observables = get(kwargs, :observables, (_store_state, ))
    G = genfunc(tlist, 1; state=state, backwards=backwards,
                storage=storage, observables=observables, init=true)
    wrk = initpropwrk(state, tlist, method, G; kwargs...)
    return propagate(state, genfunc, tlist, wrk; kwargs...)
end


function propagate(state, genfunc, tlist, method::Symbol; kwargs...)
    return propagate(state, genfunc, tlist, Val(method); kwargs...)
end


function propagate(state, genfunc, tlist, wrk;
                   backwards=false,
                   storage=nothing,
                   observables=(_store_state, ),
                   hook=nothing,
                   kwargs...)
    return_storage = false
    if storage === true
        storage = init_storage(state, tlist, observables)
        return_storage = true
    end
    state = copy(state)
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
        propstep!(state, generator, (backwards ? -dt : dt), wrk; kwargs...)
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


function propagate(state, genfunc, tlist, method::Val{:cheby}; kwargs...)
    backwards = get(kwargs, :backwards, false)
    storage = get(kwargs, :storage, nothing)
    observables = get(kwargs, :observables, (_store_state, ))
    # TODO: generate multiple examplary G
    G = genfunc(tlist, 1; state=state, backwards=backwards,
                storage=storage, observables=observables, init=true)
    wrk = initpropwrk(state, tlist, method, G; kwargs...)
    return propagate(state, genfunc, tlist, wrk; kwargs...)
end
