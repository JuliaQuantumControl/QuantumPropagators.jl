module QuantumPropagators

include("./cheby.jl")
include("./newton.jl")
include("./expprop.jl")

export cheby_coeffs, cheby_coeffs!, ChebyWrk, cheby!
export NewtonWrk, newton!
export ExpPropWrk, expprop!
export initpropwrk, init_storage, propstep!, propagate


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


"""Return a `storage` object for [`propagate`](@ref).

```julia
storage = init_storage(state, tlist, observables=(state->copy(state), ))
```

Return an Array suitable for storing the result of applying the given
`observables` to a `state` for every point in `tlist`. There are three cases:

1.  For a single observable, the resulting `storage` is a Vector of the same
    length as `tlist` that holds the result of applying the observable. For the
    default observable `state->copy(state)`, that is a Vector of states
    (`Vector{typeof(state)}`).

2.  For multiple observables, if the observables are uniform (all return an
    object of the same type), the resulting `storage` will be an Array of size
    `n × nt` where `n` is the number of `observables` and `nt` is the length of
    `tlist`.  This applies to e.g. the case where the observables are
    normal expectation values.

    ~~~julia
    observables=(state->dot(state, Ô₁, state), state->dot(state, Ô₂, state))
    ~~~

    for two Hermitian operators `Ô₁`, `Ô₂` would result in a `storage` of type
    `Matrix{Float64}`. After a propagation with [`propagate`](@ref), the
    expectation values of `Ô₁` over time would then be accessible as
    `storage[1,:]`

3.  For multiple observables, if the observables are not uniform, the resulting
    `storage` will be a Vector of length `nt` for the observable-tuples. For
    example,

    ~~~julia
    observables=(state->dot(state, Ô₁, state), state->count_poplevels(state))
    ~~~

    where `count_poplevels` is a function that counts the number of levels with
    non-zero population, the resulting `storage` would be a
    `Vector{Tuple{Float64, Int64}}`.
"""
function init_storage(state, tlist, observables=(state->copy(state), ))
    val_tuple = Tuple(O(state) for O in observables)
    first = val_tuple[1]
    nt = length(tlist)
    if length(val_tuples) == 1
        return Vector{typeof(first)}(undef, nt)
    else
        is_uniform = all(typeof(v) == typeof(first) for v in val_tuple[2:end])
        if is_uniform
            n = length(observables)
            return Array{typeof(first)}(undef, n, nt)
        else
            return Vector{typeof(val_tuple)}(undef, nt)
        end
    end
end


# helper function for `propagate`
function _prop_store(storage, i, state, observables, in_place=false)
    val_tuple = Tuple(O(state) for O in observables)
    if length(val_tuple) == 1
        if in_place
            copyto!(storage[i], val_tuple[1])
        else
            storage[i] = val_tuple[1]
        end
    else
        if length(size(storage)) == 2  # 2D Array
            for i_obs in 1:length(val_tuple)
                if in_place
                    copyto!(storage[i_obs, i], val_tuple[i_obs])
                else
                    storage[i_obs, i] = val_tuple[i_obs]
                end
            end
        else  # Vector of tuples
            if in_place
                copyto!(storage[i], val_tuple)
            else
                storage[i] = val_tuple
            end
        end
    end
end


"""Propagate a state over an entire time grid.

```julia
propagate(state, genfunc, tlist, backwards=false;
          wrk=nothing, storage=nothing, observables=(state->copy(state), ),
          storage_in_place=false)
```

propagates `state` over the time grid in `tlist`, using piecewise-constant
dynamical generators (Hamiltonians or Liouvillians) determined by `genfunc`,
and returns the resulting propagated state.

For the i'th time interval, `genfunc(tlist, i)` must return the generator for
that time interval. Generally, when approximating a time-continuous dynamical
generator as piecewise-constant on the time grid, it should be evaluated at the
*midpoint* of the interval. A possible exception is the first and last
interval, which may be better evaluated at `tlist[1]` and `tlist[end]` to
ensure exact boundary conditions like control fields that are exactly zero.

The propagation method is determined by `wrk`, see [`initpropwrk`](@ref). If
`wrk` is not given, `propagate` will attempt to choose the most appropriate
method for the given parameters.

In general, there is no requirement that `tlist` has a constant time step,
although some propagation methods (most notably [`cheby!`](@ref)) only support
a uniform time grid.

If `storage` is given as an Array, it will be filled with data determined by
the `observables`. The default "observable", `state->copy(state)` results in
the propagated states at every point in time being stored. Other use cases
would include the storage of expectation values, e.g. with

~~~julia
observables=(state->dot(state, Ô₁, state), state->dot(state, Ô₂, state))
~~~

where `Ô₁`, `Ô₂` are two Hermitian operators.

The `storage` array should be initialized with [`init_storage`](@ref). See its
documentation for the required layout of `storage` for different types of
`observables`. If `storage_in_place` is `true`, `copyto!` will be used to store
values into `storage`. This only works in special cases (like the default
storage of propagated states) and may require additional initialization beyond
[`init_storage`](@ref) (e.g., pre-allocating the states in `storage`).  Use
with caution.

Note that the term "observables" is used very loosely here: the `observables`
are not required to yield real values, but are allowed to directly
calculate, e.g., the complex amplitude α of a coherent state in quantum optics,
the number of levels with non-zero population, or the propagated state
transformed from a moving frame to a lab frame.

If `backwards` is `true`, the input state is assumed to be at time
`tlist[end]`, and the propagation progresses backwards in time (with a negative
time step `dt`). If `storage` is given, it will be filled back-to-front during
the backwards propagation.

The `propagate` routine returns the propagated state at `tlist[end]`,
respectively `tlist[1]` if `backwards=true`.
"""
function propagate(state, genfunc, tlist, backwards=false;
                   wrk=nothing,
                   storage=nothing,
                   observables=(state->copy(state), ),
                   storage_in_place=false,
                  )
    state = copy(state)
    if wrk == nothing
        wrk = initipropwrk(state, tlist, genfunc(tlist, 1))
    end
    intervals = enumerate(tlist[2:end])
    if backwards
        intervals = Iterators.reverse(intervals)
        if storage ≠ nothing
            _prop_store(storage, lastindex(storage), state, observables,
                        storage_in_place)
        end
    else
        if storage ≠ nothing
            _prop_store(storage, 1, state, observables, storage_in_place)
        end
    end
    # TODO: optional progress meter
    for (i, t_end) in intervals
        dt = t_end - tlist[i-1]
        generator = genfunc(tlist, i)
        propstep!(state, generator, (backwards ? -dt : dt), wrk)
        if storage ≠ nothing
            _prop_store(
                storage, i + (backwards ? -1 : 1), state, observables,
                storage_in_place
            )
        end
    end
    return state
end

end
