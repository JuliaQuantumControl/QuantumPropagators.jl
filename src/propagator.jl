using Printf

"""Abstract base type for all `Propagator` objects.

All `Propagator` objects must be instantiated via [`init_prop`](@ref) and
implement the following interface.

# Properties

* `state` (read-only): The current quantum state in the propagation
* `tlist` (read-only): The time grid for the propagation
* `t` (read-only): The time at which `state` is defined. An element of `tlist`.
* `parameters`: parameters that determine the dynamics. The structure of the
  parameters depends on the concrete `Propagator` type (i.e., the propagation
  method). Mutating the `parameters` affects subsequent propagation steps.
* `backward`: Boolean flag to indicate whether the propagation moves forward or
  backward in time
* `inplace`: Boolean flag to indicate whether `propagator.state` is modified
  in-place or is recreated by every call of `prop_step!` or `set_state!`.

Concrete `Propagator` types may have additional properties or fields, but these
should be considered private.

# Methods

* [`reinit_prop!`](@ref) — reset the propagator to a new initial state at the
  beginning of the time grid (or the end, for backward propagation)
* [`prop_step!`](@ref) – advance the propagator by one step forward (or
  backward) on the time grid.
* [`set_state!`](@ref) — safely mutate the current quantum `state` of the
  propagation. Note that directly mutating the `state` property is not safe.
  However, `Ψ = propagator.state; foo_mutate!(Ψ), set_state!(propagator, Ψ)`
  for some mutating function `foo_mutate!` is guaranteed to be safe and
  efficient for both in-place and not-in-place propagators.
* [`set_t!`](@ref QuantumPropagators.set_t!) — safely mutate the current time
  (`propagator.t`), snapping to the values of `tlist`.

# See also

* [`PiecewisePropagator`](@ref) — specialization of `AbstractPropagator` for
  piecewise propagation methods.
* [`PWCPropagator`](@ref) — specialization of [`PiecewisePropagator`](@ref) for
  piecewise-constant propagation methods.
"""
abstract type AbstractPropagator end

"""[`AbstractPropagator`](@ref) sub-type for piecewise propagators.

A piecewise propagator is determined by a single parameter per control and time
grid interval. Consequently, the `propagator.parameters` are a dictionary
mapping the controls found in the `generator` via [`get_controls`](@ref
QuantumPropagators.Controls.get_controls) to a vector of values defined on the
intervals of the time grid, see [`discretize_on_midpoints`](@ref). This does
not necessarily imply that these values are the piecewise-constant amplitudes
for the intervals. A general piecewise propagator might use interpolation to
obtain actual amplitudes within any given time interval.

When the amplitudes *are* piecewise-constant, the propagator should be a
concrete incantation of a [`PWCPropagator`](@ref).
"""
abstract type PiecewisePropagator <: AbstractPropagator end

"""[`PiecewisePropagator`](@ref) sub-type for piecewise-constant propagators.

Like the more general [`PiecewisePropagator`](@ref), this is characterized by
`propagator.parameters` mapping the controls in the `generator` to a vector of
amplitude value on the midpoints of the time grid intervals.

The propagation will use these values as constant within each interval.
"""
abstract type PWCPropagator <: PiecewisePropagator end


function Base.getproperty(propagator::AbstractPropagator, name::Symbol)
    # Propagators where :state, :parameters, and :t are not direct fields will
    # have to override `getproperty`
    if name == :generator
        # We do not want people to mutate an existing generator or tlist, so we
        # actively present access to these fields
        error("type $(nameof(typeof(propagator))) has no property $(name)")
    end
    return getfield(propagator, name)
end

function Base.setproperty!(propagator::AbstractPropagator, name::Symbol, value)
    if name ≡ :state
        # A user would probably expect the following behavior:
        #
        #     propagator.state = state;
        #     @assert propagator.state === state
        #
        # At least for in-place propagators, this expectation will not hold,
        # since we'll want to use `copyto!` to overwrite the internal state of
        # the `propagator` instead of changing the reference. To make the
        # behaviour explicit, we for the user to use `set_state!`
        error("The state of a propagator can only be set via `set_state!`")
    elseif name ≡ :parameters
        # Optimal control is most straightforward to implement if
        # propagator.parameters can alias an external (mutable) data structure.
        # Thus, it is best to give direct access to this property
        setfield!(propagator, :parameters, value)
    elseif name ≡ :tlist
        error("The tlist of a propagator is read-only")
    elseif name ≡ :t
        # The time must always be an element of `tlist`, so setting it to an
        # arbitrary value is not possible
        error("The curerent time of a propagator can only be set via `set_t!`")
    elseif name ≡ :generator
        # See `getproperty`
        error("type $(nameof(typeof(propagator))) has no property $(name)")
    else
        setfield!(propagator, name, value)
    end
end

function Base.propertynames(propagator::AbstractPropagator, private::Bool=false)
    public_properties = (:state, :tlist, :t, :parameters, :backward, :inplace)
    if private
        return Tuple(union(public_properties, fieldnames(typeof(propagator))))
    else
        return public_properties
    end
end


"""Initialize a `Propagator`.

```julia
propagator = init_prop(
    state, generator, tlist;
    method=:auto,
    backward=false,
    inplace=true,
    piecewise=nothing,
    pwc=nothing,
    kwargs...
)
```

initializes a propagator for the time propagation of the given `state` over a
time grid `tlist` under the time-dependent generator (Hamiltonian/Liouvillian)
`generator`.

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
  given `state`, `tlist`, and `generator`, cf. [`choose_propmethod`](@ref)
* `backward`: If `true`, initialize the propagator for a backward propagation.
  The resulting `propagator.t` will be `tlist[end]`, and subsequent calls to
  [`prop_step!`](@ref) will move backward on `tlist`.
* `inplace`: If `true`, the `state` property of the resulting propagator will
  be changed in-place by any call to [`prop_step!`](@ref). If `false`, each
  call to [`prop_step!`](@ref) changes the reference for `propagator.state`,
  and the propagation will not use any in-place operations. Not all propagation
  methods may support both in-place and not-in-place propagation. In-place
  propagation is generally more efficient but may not be compatible, e.g., with
  automatic differentiation.
* `piecewise`: If given as a boolean, `true` enforces that the resulting
  propagator is a [`PiecewisePropagator`](@ref), and `false` enforces is not to
  be a [`PiecewisePropagator`](@ref)
* `pwc`: Like `piecewise`, for for the stronger [`PWCPropagator`](@ref)

All other `kwargs` are method-dependent and are ignored for methods that do not
support them.

The type of the returned `propagator` is a sub-type of
[`AbstractPropagator`](@ref), respectively a sub-type of
[`PiecewisePropagator`](@ref) if `piecewise=true` or a sub-type of
[`PWCPropagator`](@ref) if `pwc=true`.

# See also

* [`reinit_prop!`](@ref) — Re-initialize a propagator
* [`propagate`](@ref) — Higher-level propagation interface
* `QuantumPropagators.Interfaces.check_propagator` — a function to verify the
  interface described above.
"""
function init_prop(
    state,
    generator,
    tlist,
    method::Val{:auto};
    backward=false,
    inplace=true,
    verbose=false,
    piecewise=nothing,
    pwc=nothing,
    kwargs...
)
    method = choose_propmethod(generator, state, tlist; piecewise, pwc, inplace)
    if verbose
        @info "Auto-choosing propagation method $(repr(typeof(method).parameters[1]))"
    end
    propagator = init_prop(state, generator, tlist, method; backward=backward, kwargs...)
    if piecewise ≡ true
        if !(propagator isa PiecewisePropagator)
            error("Cannot initalize propagator as piecewise with method=$(repr(method))")
        end
    end
    if pwc ≡ true
        if !(propagator isa PWCPropagator)
            error(
                "Cannot initalize propagator as piecewise-constant with method=$(repr(method))"
            )
        end
    end
    return propagator
end

function init_prop(state, generator, tlist, method; kwargs...)
    throw(
        TypeError(
            :init_prop,
            "`method` must be a known symbol, not $(repr(method))",
            Symbol,
            typeof(method)
        )
    )
end


# `method` as a keyword argument
function init_prop(state, generator, tlist; method=Val(:auto), kwargs...)
    return init_prop(state, generator, tlist, method; kwargs...)
end


# `method` as a symbol
function init_prop(state, generator, tlist, method::Symbol; kwargs...)
    return init_prop(state, generator, tlist, Val(method); kwargs...)
end


function _get_uniform_dt(tlist::Vector; tol=1e-12, warn=false)
    dt = float(tlist[2] - tlist[1])
    for i = 2:length(tlist)-1
        dt_i = tlist[i+1] - tlist[i]
        Δ = abs(dt_i - dt)
        if Δ > tol
            if warn
                @warn "Non-uniform time grid: dt = $(@sprintf("%.2e", dt_i)) in interval $i differs from the first dt=$(@sprintf("%.2e", dt)) by Δ = $(@sprintf("%.2e", Δ)) > tol = $(@sprintf("%.2e", tol))"
            end
            return nothing
        end
    end
    return dt
end


"""Choose a suitable propagation method.

```julia
method = choose_propmethod(generator, state, tlist;
                           pwc=nothing, piecewise=nothing, inplace=true)
```

identifies a suitable propagation method for the given `generator`, `state` and
`tlist`. If `piecewise` or `pwc` are given as `true`, only consider methods
that result in in a [`PiecewisePropagator`](@ref) or [`PWCPropagator`](@ref),
respectively. If `piecewise` or `pwc` are given as `false`, disregard any
methods that result in these propagators. Only propagators that support the
given `inplace` are taken into account.
"""
function choose_propmethod(
    generator,
    state,
    tlist;
    pwc=nothing,
    piecewise=nothing,
    inplace=true
)
    (pwc ≡ false) && error("No non-pwc propagator available")
    (piecewise ≡ false) && error("No non-piecewise propagator available")
    method = nothing
    is_small = false
    try
        is_small = length(state) <= 8
    catch exception
        if isa(exception, MethodError)
            @warn "Cannot determine problem dimension"
        else
            rethrow()
        end
    end
    if ((pwc ≢ false) || (piecewise ≢ false))
        if is_small
            method = Val(:expprop)
        else
            if inplace
                method = Val(:newton)
            end
            controls = get_controls(generator)
            genop = _pwc_get_max_genop(generator, controls, tlist)
            if has_real_eigvals(genop) ≡ true
                dt = _get_uniform_dt(tlist)
                if !isnothing(dt)
                    method = Val(:cheby)
                end
            end
        end
    end
    if isnothing(method)
        error("Cannot find a suitable propagation method.")
    end
    return method
end


"""Re-initialize a propagator.

```julia
reinit_prop!(propagator, state; kwargs...)
```

resets the `propagator` to `state` at the beginning of the time grid,
respectively the end of the time grid if `propagator.backward` is true.

At a minimum, this is equivalent to a call to [`set_state!`](@ref) follow by a
call to [`set_t!`](@ref), but some propagators may have additional requirements
on re-initialization, such as refreshing expansion coefficients for
[`ChebyPropagator`](@ref). In this case, the `kwargs` may be additional keyword
arguments specific to the concrete type of propagator.
"""
function reinit_prop!(propagator, state; kwargs...)
    set_state!(propagator, state)
    if propagator.backward
        set_t!(propagator, propagator.tlist[end])
    else
        set_t!(propagator, propagator.tlist[begin])
    end
end


"""Advance the `propagator` by a single time step.

```julia
state = prop_step!(propagator)
```

returns the state obtained from propagating to the next point on the time grid
from `propagator.t`, respectively the previous point if `propagator.backward`
is true.

When the propagation would lead out of the time grid, `prop_step!` leaves
`propagator` unchanged and returns `nothing`. Thus, a return value of `nothing`
may be used to signal that a propagation has completed.
"""
function prop_step! end

"""Set the current `state` of the `propagator`.

```julia
set_state!(propagator, state)
```

sets the `propagator.state` property. In order to mutate the current state
after a call to [`prop_step!`](@ref), the following pattern is recommended:

```
Ψ = propagator.state
foo_mutate!(Ψ)
set_state!(propagator, Ψ)
```

where `foo_mutate!` is some function that mutates `Ψ`.  This is guaranteed to
work efficiently both for in-place and not-in-place propagators, without
incurring unnecessary copies.

!!! warning
    ```
    foo_mutate!(propagator.state)
    ```

    by itself is not a safe operation. Always follow it by

    ```
    set_state!(propagator, propagator.state)
    ```

# See also

* [`set_t!`](@ref QuantumPropagators.set_t!) — set `propagator.t`.

"""
function set_state!(propagator::AbstractPropagator, state)
    if state ≢ propagator.state
        if propagator.inplace
            copyto!(propagator.state, state)
        else
            setfield!(propagator, :state, state)
        end
    end
end


"""Set the current time for the propagation.

```julia
set_t!(propagator, t)
```

Sets `propagator.t` to the given value of `t`, where `t` must be an element of
`propagator.tlist`.

# See also

* [`set_state!`](@ref) — set `propagator.state`
"""
function set_t!(propagator, t)
    tlist = propagator.tlist
    if t <= tlist[1]
        n = 1
    else
        N = length(tlist)
        if t >= tlist[end]
            n = N
        else
            n = min(searchsortedfirst(tlist, t), N)
        end
    end
    (t ≈ tlist[n]) || (@warn ("Snapping t=$t to time grid value $(tlist[n])"))
    setfield!(propagator, :t, tlist[n])
end
