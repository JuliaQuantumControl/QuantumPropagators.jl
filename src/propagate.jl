using LinearAlgebra
using ProgressMeter

using .Storage: init_storage, write_to_storage!
import .Storage: map_observables

using .Interfaces: check_state, check_generator


# Return `true` if `H` has only real eigenvalues, `false` if `H` has
# eigenvalues with a non-zero imaginary part, and `nothing` if field of the
# eigenvalues of `H` cannot be determined
function has_real_eigvals(H)
    if hasmethod(ishermitian, (typeof(H),))
        return ishermitian(H)
    else
        @warn "ishermitian is not defined for $(typeof(H))"
        return nothing
    end
end


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
    method=:auto,
    check=true,
    backward=false,
    inplace=true,
    verbose=false,
    piecewise=nothing,
    pwc=nothing,
    storage=nothing,
    observables=<store state>,
    callback=nothing,
    showprogress=false,
    init_prop_kwargs...)
```

propagates `state` of the entire time grid and returns the propagates states,
or a storage array of data collected during the propagation.

# Arguments

* `state`: The "initial" state for the propagation. For `backward=false`, this
  state is taken to be at initial time (`tlist[begin]`); and for
  `backward=true`, at the final time (`tlist[end]`)
* `generator`: The time-dependent generator of the dynamics
* `tlist`: The time grid over which which the propagation is defined. This may
  or may not be equidistant.

# Keyword arguments

* `method`: The propagation method to use. The default value of `:auto`
  attempts to choose the best method available, based on the properties of the
  given `state`, `tlist`, and `generator`.
* `check`: if `true`, check that `state` and `generator` pass
  [`check_state`](@ref) and [`check_generator`](@ref).
* `backward`: If `true`, propagate backward in time
* `inplace`: If `true`, propagate using in-place operations. If `false`, avoid
  in-place operations. Not all propagation methods
  support both in-place and not-in-place propagation.
* `piecewise`: If given a a boolean, limit the propagation to "piecewise"
  methods, respectively disallow piecewise methods
* `pwc`: If given a a boolean, limit the propagation to piecewise-constant
  methods, respectively disallow piecewise-constant methods
* `storage`: Flag whether to store and return the propagated states /
  observables, or pre-allocated storage array. See Notes below.
* `observables`: Converters for data to be stored in `storage`. See Notes
  below.
* `callback`: Function to call after each propagation step. See Notes below.
* `showprogess`: Whether to show a progress bar. See Notes below.

All remaining keyword arguments are passed to [`init_prop`](@ref) to initialize
the [`Propagator`](@ref AbstractPropagator) that is used internally to drive
the optimization. Unknown keyword arguments will be ignored.

# Notes

In general, there is no requirement that `tlist` has a constant time step,
although some propagation methods (most notably [`cheby!`](@ref
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
    state, generator, tlist; method=:auto
    backward=false; storage=true, observables=observables,
    callback=nothing, showprogress=false, init_prop_kwargs...)
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

If `showprogress` is given as `true`, a progress bar will be shown for
long-running propagationn. In order to customize the progress bar,
`showprogress` may also be a function that receives `length(tlist)` and returns
a `ProgressMeter.Progress` instance.

If `in_place=false` is given, the propagation avoids in-place operations. This
is slower than `inplace=true`, but is often required in the context of
automatic differentiation (AD), e.g., with
[Zygote](https://fluxml.ai/Zygote.jl/). That is, use `in_place=false` if
`propagate` is called inside a function to be passed to `Zygote.gradient`,
`Zygote.pullback`, or a similar function. In an AD context, `storage` and
`showprogress` should not be used.

The `propagate` routine returns the propagated state at `tlist[end]`,
respectively `tlist[1]` if `backward=true`, or a storage array with the
stored states / observable data if `storage=true`.

# See also

* [`init_prop`](@ref) — Propagate via a [`Propagator`](@ref AbstractPropagator)
  object
"""
function propagate(state, generator, tlist; method=Val(:auto), kwargs...)
    return propagate(state, generator, tlist, method; kwargs...)
end


function propagate(state, generator, tlist, method::Symbol; kwargs...)
    return propagate(state, generator, tlist, Val(method); kwargs...)
end


# default `propagate` implementation
function propagate(
    state,
    generator,
    tlist,
    method::Val;
    check=true,
    storage=nothing,
    showprogress=false,
    observables=_StoreState(),
    callback=nothing,
    inplace=true, # cf. default of init_prop
    for_expval=true,  # undocumented
    for_immutable_state=true,  # undocumented
    for_mutable_state=inplace,  # undocumented
    kwargs...
)
    backward = get(kwargs, :backward, false)
    atol = get(kwargs, :atol, 1e-14)  # for checks
    quiet = get(kwargs, :quiet, false)  # for checks

    if check
        if !(tlist isa Vector{Float64})
            error(
                "The `tlist` in `propagate` must be a Vector{Float64}, not $(typeof(tlist))"
            )
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

    propagator = init_prop(state, generator, tlist, method; inplace, kwargs...)

    return_storage = false
    if storage === true
        storage = init_storage(state, tlist, observables)
        return_storage = true
    end
    state = copy(state)
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
    if isa(showprogress, Function)
        progressmeter = showprogress(N)
    else
        progressmeter = Progress(N, enabled=showprogress)
        if !showprogress
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
