using LinearAlgebra
using ProgressMeter

using .Storage: init_storage, write_to_storage!
import .Storage: map_observables

using .Interfaces: check_state, check_generator, check_tlist


# Default "observables" for just storing the propagated state. We want to keep
# this private: It doesn't combine with other observables and depends on
# internals of the default init_storage.
struct _StoreState end
map_observables(::_StoreState, tlist, i, state) = copy(state)
map_observables(::_StoreState, tlist, i, state::Vector) = state


# Work around https://github.com/timholy/ProgressMeter.jl/issues/214
import ProgressMeter
struct NoProgress end
ProgressMeter.next!(p::NoProgress) = nothing


"""Propagate a state over an entire time grid.

```julia
state = propagate(
    state,
    generator,
    tlist;
    method,  # mandatory keyword argument
    check=true,
    backward=false,
    inplace=true,
    verbose=false,
    piecewise=nothing,
    pwc=nothing,
    storage=nothing,
    observables=<store state>,
    callback=nothing,
    show_progress=false,
    init_prop_kwargs...)
```

propagates `state` of the entire time grid and returns the propagated states,
or a storage array of data collected during the propagation. This high-level
routine performs the following three steps:

1.  If `check=true` (default), check that `state`, `generator`, and `tlist` are
    consistent with the required interface.

2.  Initialize a `propagator` via [`init_prop`](@ref):

    ```
    init_prop(state, generator, tlist; method, inplace, init_prop_kwargs...)
    ```

3.  Call and return the result of

    ```
    propagate(propagator; storage, observables, show_progress, callback)
    ```

# Arguments

* `state`: The "initial" state for the propagation. For `backward=false`, this
  state is taken to be at initial time (`tlist[begin]`); and for
  `backward=true`, at the final time (`tlist[end]`)
* `generator`: The time-dependent generator of the dynamics
* `tlist`: The time grid over which which the propagation is defined. This may
  or may not be equidistant.

# Mandatory keyword arguments

* `method`: The propagation method to use. May be given as a name
  (`Symbol`), but the recommended usage is to pass a module implementing the
  propagation method, cf. [`init_prop`](@ref).

# Optional keyword arguments

* `check`: if `true`, check that `state`, `generator`, and `tlist` pass
  [`check_state`](@ref), [`check_generator`](@ref) and [`check_tlist`](@ref),
  respectively.
* `backward`: If `true`, propagate backward in time
* `inplace`: If `true`, propagate using in-place operations. If `false`, avoid
  in-place operations. Not all propagation methods
  support both in-place and not-in-place propagation.
* `piecewise`: If given as a boolean, ensure that the internal `propagator` is
  an instance of [`PiecewisePropagator`](@ref), cf. [`init_prop`](@ref).
* `pwc`: If given a a boolean, do a piecewise constant propagation where the
  generator in each interval is constant (the internal `propagator` is a
  [`PWCPropagator`](@ref), cf. [`init_prop`](@ref))
* `storage`: Flag whether to store and return the propagated states /
  observables, or pre-allocated storage array. See Notes below.
* `observables`: Converters for data to be stored in `storage`. See Notes
  below.
* `callback`: Function to call after each propagation step. See Notes below.
* `show_progress`: Whether to show a progress bar. See Notes below.

All remaining keyword arguments are passed to [`init_prop`](@ref) to initialize
the [`Propagator`](@ref AbstractPropagator) that is used internally to drive
the optimization. Unknown keyword arguments will be ignored.

# Notes

In general, there is no requirement that `tlist` has a constant time step,
although some propagation methods (most notably [`Cheby`](@ref
QuantumPropagators.Cheby.cheby!)) only support a uniform time grid.

If `storage` is given as a container pre-allocated via [`init_storage`](@ref),
it will be filled with data determined by the `observables`. Specifically,
after each propagation step,

```julia
data = map_observables(observables, tlist, i, state)
write_to_storage!(storage, i, data)
```

is executed, where `state` is defined at time `tlist[i]`.
See [`map_observables`](@ref) and [`write_to_storage!`](@ref) for details. The
default values for `observables` results simply in the propagated states at
every point in time being stored.

The `storage` parameter may also be given as `true`, and a new storage array
will be created internally with [`init_storage`](@ref) and returned instead of
the propagated state:

```julia
data = propagate(
    state, generator, tlist; method,
    backward=false; storage=true, observables=observables,
    callback=nothing, show_progress=false, init_prop_kwargs...)
```

If `backward` is `true`, the input state is assumed to be at time
`tlist[end]`, and the propagation progresses backward in time (with a negative
time step `dt`). If `storage` is given, it will be filled back-to-front during
the backward propagation.

If `callback` is given as a callable, it will be called after each propagation
step, as `callback(propagator, observables)` where `propagator` is
[`Propagator`](@ref AbstractPropagator) object driving the propagation.
The `callback` is called before calculating any observables.
Example usage includes writing data to file, or modifying `state` via
[`set_state!`](@ref), e.g., removing amplitude from the lowest and highest
level to mitigate "truncation error".

If `show_progress` is given as `true`, a progress bar will be shown for
long-running propagation. In order to customize the progress bar,
`show_progress` may also be a function that receives `length(tlist)` and
returns a `ProgressMeter.Progress` instance.

If `in_place=false` is given, the propagation avoids in-place operations. This
is slower than `inplace=true`, but is often required in the context of
automatic differentiation (AD), e.g., with
[Zygote](https://fluxml.ai/Zygote.jl/). That is, use `in_place=false` if
`propagate` is called inside a function to be passed to `Zygote.gradient`,
`Zygote.pullback`, or a similar function. In an AD context, `storage` and
`show_progress` should not be used.

The `propagate` routine returns the propagated state at `tlist[end]`,
respectively `tlist[1]` if `backward=true`, or a storage array with the
stored states / observable data if `storage=true`.
"""
function propagate(
    state,
    generator,
    tlist;
    method,
    check=true,
    storage=nothing,
    show_progress=false,
    observables=_StoreState(),
    callback=nothing,
    inplace=true, # cf. default of init_prop
    for_expval=true,  # undocumented
    for_immutable_state=true,  # undocumented
    for_mutable_operator=inplace,  # undocumented
    for_mutable_state=inplace,  # undocumented
    kwargs...
)
    atol = get(kwargs, :atol, 1e-14)  # for checks
    quiet = get(kwargs, :quiet, false)  # for checks

    if check
        if !(tlist isa Vector{Float64})
            error(
                "The `tlist` in `propagate` must be a Vector{Float64}, not $(typeof(tlist))"
            )
        end
        valid_tlist =
            check_tlist(tlist; quiet, _message_prefix="On `tlist` in `propagate`: ")
        if !valid_tlist
            error("The `tlist` in `propagate` does not pass `check_tlist`")
        end
        valid_state = check_state(
            state;
            for_immutable_state,
            for_mutable_state,
            atol,
            quiet,
            _message_prefix="On initial `state` in `propagate`: "
        )
        if !valid_state
            error("The initial state in `propagate` does not pass `check_state`")
        end
        valid_generator = check_generator(
            generator;
            state=state,
            tlist=tlist,
            for_immutable_state,
            for_mutable_operator,
            for_mutable_state,
            for_expval,
            atol,
            quiet,
            _message_prefix="On `generator` in `propagate`: "
        )
        if !valid_generator
            error(
                "The generator $(repr(generator)) in `propagate` does not pass `check_generator`"
            )
        end
    end

    propagator = init_prop(state, generator, tlist; method, inplace, kwargs...)
    # Calling `init_prop` with `method` as a keyword argument instead of a
    # positional arguments allows for overriding `init_prop` for some type of
    # `state` and/or `generator`, independent of the `method`.
    # See `test_nonstandard_generators.jl`

    return propagate(propagator; storage, observables, show_progress, callback)

end


"""
```julia
state = propagate(
    state,
    propagator;
    storage=nothing,
    observables=<store state>,
    show_progress=false,
    callback=nothing,
    reinit_prop_kwargs...
)
```

re-initializes the given `propagator` with `state` (see [`reinit_prop!`](@ref))
and then calls the lower-level `propagate(propagator; ...)`.
"""
function propagate(
    state,
    propagator;
    storage=nothing,
    observables=_StoreState(),
    show_progress=false,
    callback=nothing,
    kwargs...
)
    reinit_prop!(propagator, state; kwargs...)
    return propagate(propagator; storage, observables, show_progress, callback)
end


"""
```julia
state = propagate(
    propagator;
    storage=nothing,
    observables=<store state>,
    show_progress=false,
    callback=nothing,
)
```

propagates a freshly initialized `propagator` (immediately after
[`init_prop`](@ref)). Used in the higher-level
[`propagate(state, generator, tlist; kwargs...)`](@ref).
"""
function propagate(
    propagator;
    storage=nothing,
    observables=_StoreState(),
    show_progress=false,
    callback=nothing,
)

    state = propagator.state
    tlist = propagator.tlist
    backward = propagator.backward
    return_storage = false
    if storage === true
        storage = init_storage(state, tlist, observables)
        return_storage = true
    end
    intervals = enumerate(tlist[2:end])
    if backward
        intervals = Iterators.reverse(intervals)
        if storage ≠ nothing
            _write_to_storage!(storage, lastindex(tlist), state, observables, tlist)
        end
    else
        if storage ≠ nothing
            _write_to_storage!(storage, 1, state, observables, tlist)
        end
    end

    N = length(intervals)
    if isa(show_progress, Function)
        progressmeter = show_progress(N)
    else
        progressmeter = Progress(N, enabled=show_progress)
        if !show_progress
            # XXX: https://github.com/timholy/ProgressMeter.jl/issues/214
            progressmeter = NoProgress()
        end
    end

    for (i, t_end) in intervals
        prop_step!(propagator)
        if storage ≠ nothing
            _write_to_storage!(
                storage,
                i + (backward ? 0 : 1),
                propagator.state,
                observables,
                tlist
            )
        end
        if callback ≠ nothing
            callback(propagator, observables)
        end
        next!(progressmeter)
    end
    if return_storage
        return storage
    else
        return propagator.state
    end

end

@inline function _write_to_storage!(storage, i::Integer, state, observables, tlist)
    # This specific implementation is referenced in the documentation of
    # `propagate`. It should not be customized.
    data = map_observables(observables, tlist, i, state)
    write_to_storage!(storage, i, data)
end
