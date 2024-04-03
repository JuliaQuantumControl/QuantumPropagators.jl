module Controls

export discretize, discretize_on_midpoints
export get_controls, get_parameters
export ParameterizedFunction
export get_tlist_midpoints
export evaluate, evaluate!, substitute

using LinearAlgebra


"""Evaluate `control` at every point of `tlist`.

```julia
values = discretize(control, tlist; via_midpoints=true)
```

discretizes the given `control` to a Vector of values defined on the points of
`tlist`.

If `control` is a function, it is first evaluated at the midpoint of `tlist`,
see [`discretize_on_midpoints`](@ref), and then the values on the midpoints are
converted to values on `tlist`. This discretization is more stable than
directly evaluating the control function at the values of `tlist`, and ensures
that repeated round-trips between [`discretize`](@ref) and
[`discretize_on_midpoints`](@ref) can be done safely, see the note in the
documentation of [`discretize_on_midpoints`](@ref).

The latter can still be achieved by passing `via_midpoints=false`. While such a
direct discretization is suitable e.g. for plotting, but it is unsuitable
for round-trips between [`discretize`](@ref) and
[`discretize_on_midpoints`](@ref)  (constant controls on `tlist` may result in
a zig-zag on the intervals of `tlist`).

If `control` is a vector, a copy of `control` will be returned if it is of the
same length as `tlist`. Otherwise, `control` must have one less value than
`tlist`, and is assumed to be defined on the midpoints of `tlist`. In that
case, [`discretize`](@ref) acts as the inverse of
[`discretize_on_midpoints`](@ref). See [`discretize_on_midpoints`](@ref) for
how control values on `tlist` and control values on the intervals of `tlist`
are related.
"""
function discretize(control::Function, tlist; via_midpoints=true)
    if via_midpoints
        vals_on_midpoints = discretize_on_midpoints(control, tlist)
        return discretize(vals_on_midpoints, tlist)
    else
        return [control(t) for t in tlist]
    end
end

function discretize(control::Vector, tlist)
    if length(control) == length(tlist)
        return copy(control)
    elseif length(control) == length(tlist) - 1
        # convert `control` on intervals to values on `tlist`
        # cf. pulse_onto_tlist in Python krotov package
        vals = zeros(eltype(control), length(control) + 1)
        vals[1] = control[1]
        vals[end] = control[end]
        for i = 2:length(vals)-1
            vals[i] = 0.5 * (control[i-1] + control[i])
        end
        return vals
    else
        throw(ArgumentError("control array must be defined on intervals of tlist"))
    end
end


"""Shift time grid values to the interval midpoints

```julia
tlist_midpoints = get_tlist_midpoints(
    tlist; preserve_start=true, preserve_end=true
)
```

takes a vector `tlist` of length ``n`` and returns a `Vector{Float64}` of
length ``n-1`` containing the midpoint values of each interval. The intervals
in `tlist` are not required to be uniform.

By default, the first and last point of `tlist` is preserved, see
[`discretize_on_midpoints`](@ref). This behavior can be disabled by passing
`preserve_start` and `preserve_end` as `false` in order to use the midpoints of
the first and last interval, respectively.
"""
function get_tlist_midpoints(tlist::AbstractVector; preserve_start=true, preserve_end=true)
    N = length(tlist)
    if N < 3
        @error "N = length(tlist) < 3" N
        msg = "In `get_tlist_midpoints`, argument `tlist` must have a length of at least 3"
        throw(ArgumentError(msg))
    end
    tlist_midpoints = zeros(Float64, N - 1)
    if preserve_start
        tlist_midpoints[1] = tlist[begin]
    else
        dt = float(tlist[begin+1] - tlist[begin])
        @assert dt > 0.0
        tlist_midpoints[1] = tlist[begin] + 0.5 * dt
    end
    if preserve_end
        tlist_midpoints[end] = tlist[end]
    else
        dt = float(tlist[end] - tlist[end-1])
        @assert dt > 0.0
        tlist_midpoints[end] = tlist[end-1] + 0.5 * dt
    end
    for i = 2:length(tlist_midpoints)-1
        dt = float(tlist[begin+i] - tlist[begin+i-1])
        @assert dt > 0.0
        tlist_midpoints[i] = tlist[begin+i-1] + 0.5 * dt
    end
    return tlist_midpoints
end


@doc raw"""
Evaluate `control` at the midpoints of `tlist`.

```julia
values = discretize_on_midpoints(control, tlist)
```

discretizes the given `control` to a Vector of values on the midpoints of
`tlist`. Hence, the resulting `values` will contain one less value than
`tlist`.

If `control` is a vector of values defined on `tlist` (i.e., of the same length
as `tlist`), it will be converted to a vector of values on the intervals of
`tlist`. The value for the first and last "midpoint" will remain the original
values at the beginning and end of `tlist`, in order to ensure exact boundary
conditions. For all other midpoints, the value for that midpoint will be
calculated by "un-averaging".

For example, for a `control` and `tlist` of length 5, consider the following
diagram:

~~~
tlist index:       1   2   3   4   5
tlist:             ⋅   ⋅   ⋅   ⋅   ⋅   input values cᵢ (i ∈ 1..5)
                   |̂/ ̄ ̄ ̂\ / ̂\ / ̂ ̄ ̄\|̂
midpoints:         x     x   x     x   output values pᵢ (i ∈ 1..4)
midpoints index:   1     2   3     4
~~~

We will have ``p₁=c₁`` for the first value, ``p₄=c₅`` for the last value. For
all other points, the control values ``cᵢ = \frac{p_{i-1} + p_{i}}{2}`` are the
average of the values on the midpoints. This implies the "un-averaging" for the
midpoint values ``pᵢ = 2 c_{i} - p_{i-1}``.

!!! note

    An arbitrary input `control` array may not be compatible with the above
    averaging formula. In this case, the conversion will be "lossy"
    ([`discretize`](@ref) will not recover the original `control` array; the
    difference should be considered a "discretization error"). However, any
    *further* round-trip conversions between points and intervals are bijective
    and preserve the boundary conditions. In this case, the
    [`discretize_on_midpoints`](@ref) and [`discretize`](@ref) methods are each
    other's inverse. This also implies that for an optimal control procedure,
    it is safe to modify *midpoint* values. Modifying the the values on the
    time grid directly on the other hand may accumulate discretization errors.

If `control` is a vector of one less length than `tlist`, a copy of `control`
will be returned, under the assumption that the input is already properly
discretized.

If `control` is a function, the function will be directly evaluated at the
midpoints marked as `x` in the above diagram..
"""
function discretize_on_midpoints(control::T, tlist) where {T<:Function}
    tlist_midpoints = get_tlist_midpoints(tlist)
    return discretize(control, tlist_midpoints; via_midpoints=false)
end

function discretize_on_midpoints(control::Vector, tlist)
    if length(control) == length(tlist) - 1
        return copy(control)
    elseif length(control) == length(tlist)
        vals = Vector{eltype(control)}(undef, length(tlist) - 1)
        vals[1] = control[1]
        vals[end] = control[end]
        for i = 2:length(vals)-1
            vals[i] = 2 * control[i] - vals[i-1]
        end
        return vals
    else
        throw(ArgumentError("control array must be defined on the points of tlist"))
    end
end


"""Extract a Tuple of controls.

```julia
controls = get_controls(generator)
```

extracts the controls from a single dynamical generator.

For example, if `generator = hamiltonian(H0, (H1, ϵ1), (H2, ϵ2))`, extracts
`(ϵ1, ϵ2)`.
"""
get_controls(ampl::Function) = (ampl,)
get_controls(ampl::Vector) = (ampl,)
get_controls(ampl::Number) = Tuple([])


"""
```julia
get_controls(operator)
```

for a static operator (matrix) returns an empty tuple.
"""
function get_controls(operator::AbstractMatrix)
    return Tuple([])
end



"""Evaluate all controls.

In general, `evaluate(object, args...; vals_dict=IdDict())` evaluates the
`object` for a specific point in time indicated by the positional `args`. Any
control in `object` is evaluated at the specified point in time. Alternatively,
the `vals_dict` maps a controls to value ("plug in this value for the given
control")

For example,

```julia
op = evaluate(generator, t)
```

evaluates `generator` at time `t`. This requires that any control in
`generator` is a callable that takes `t` as a single argument.

```julia
op = evaluate(generator, tlist, n)
```

evaluates `generator` for the n'th interval of `tlist`. This uses the
definitions for the midpoints in [`discretize_on_midpoints`](@ref).
The controls in `generator` may be vectors (see [`discretize`](@ref),
[`discretize_on_midpoints`](@ref)) or callables of `t`.

```julia
op = evaluate(generator, t; vals_dict)
op = evaluate(generator, tlist, n; vals_dict)
```

resolves any explicit time dependencies in `generator` at the specified point
in time, but uses the value in the given `vals_dict` for any control in
`vals_dict`.


```julia
a = evaluate(ampl, tlist, n; vals_dict=IdDict())
a = evaluate(ampl, t; vals_dict=IdDict())
```

evaluates a control amplitude to a scalar by evaluating any explicit time
dependency, and by replacing each control with the corresponding value in
`vals_dict`.

Calling `evaluate` for an object with no implicit or explicit time dependence
should return the object unchanged.

For generators without any explicit time dependence,

```julia
op = evaluate(generator; vals_dict)
```

can be used. The `vals_dict` in this case must contain values for all controls
in `generator`.

# See also:

* [`evaluate!`](@ref) — update an existing operator with a re-evaluation of a
generator at a different point in time.
"""
function evaluate(object, args...; vals_dict=IdDict())
    # Fallback: If `object` has to components, just look up `object` in
    # `vals_dict`
    return get(vals_dict, object, object)
end


function evaluate(operator::AbstractMatrix, args...; kwargs...)
    # objects without any explicit or explicit time dependency evaluate to
    # themselves
    return operator
end


# Midpoint of n'th interval of tlist, but snap to beginning/end (that's
# because any S(t) is likely exactly zero at the beginning and end, and we
# want to use that value for the first and last time interval)
function _t(tlist, n)
    @assert 1 <= n <= (length(tlist) - 1)  # n is an *interval* of `tlist`
    if n == 1
        t = tlist[begin]
    elseif n == length(tlist) - 1
        t = tlist[end]
    else
        dt = tlist[n+1] - tlist[n]
        t = tlist[n] + dt / 2
    end
    return t
end


function evaluate(func::Function, tlist::Vector, n::Int64; vals_dict=IdDict())
    if haskey(vals_dict, func)
        return vals_dict[func]
    else
        return func(_t(tlist, n))
    end
end


function evaluate(func::Function, t::Float64; vals_dict=IdDict())
    if haskey(vals_dict, func)
        return vals_dict[func]
    else
        return func(t)
    end
end


function evaluate(control::Vector, tlist::Vector, n::Int64; vals_dict=IdDict())
    if haskey(vals_dict, control)
        return vals_dict[control]
    else
        if length(control) == length(tlist) - 1
            return control[n]
        elseif length(control) == length(tlist)
            # convert to midpoint values
            if n == 1
                return control[1]
            elseif n == length(tlist)
                return control[n]
            else
                return 2 * control[n] - control[n-1]
            end
        else
            error(
                "control (length $(length(control))) must be discretized either on `tlist` (length $(length(tlist))) or on the midpoints of `tlist`"
            )
        end
    end
end


function evaluate(control::Vector, t::Float64; vals_dict=IdDict())
    if haskey(vals_dict, control)
        return vals_dict[control]
    else
        error(
            "`evaluate(control::Vector, t::Float64)` is invalid. Use e.g. `evaluate(…, tlist, n)`."
        )
    end
end


function evaluate(generator::Tuple, args...; vals_dict=IdDict())
    if isa(generator[1], Tuple)
        control = generator[1][2]
        coeff = evaluate(control, args...; vals_dict)
        if coeff isa Number
            op = coeff * generator[1][1]
        else
            error(
                "In `evaluate(::$(typeof(generator)), …)`, the control $control does not evaluate to a number with the given arguments $args, vals_dict=$vals_dict"
            )
        end
    else
        op = copy(generator[1])
    end
    for part in generator[2:end]
        if isa(part, Tuple)
            control = part[2]
            coeff = evaluate(control, args...; vals_dict)
            if coeff isa Number
                op += coeff * part[1]
            else
                error(
                    "In `evaluate(::$(typeof(generator)), …)`, the control $control does not evaluate to a number with the given arguments $args, vals_dict=$vals_dict"
                )
            end
        else
            op += part
        end
    end
    return op
end


"""Update an existing evaluation of a `generator`.

```julia
evaluate!(op, generator, args..; vals_dict=IdDict())
```

performs an in-place update on an `op` the was obtained from a previous call to
[`evaluate`](@ref) with the same `generator`, but for a different point in time
and/or different values in `vals_dict`.
"""
function evaluate!(op, generator::Tuple, args...; vals_dict=IdDict())
    if generator[1] isa Tuple
        control = generator[1][2]
        copyto!(op, generator[1][1])
        coeff = evaluate(control, args...; vals_dict)
        @assert coeff isa Number
        lmul!(coeff, op)
    else
        copyto!(op, generator[1])
    end
    for part in generator[2:end]
        if part isa Tuple
            control = part[2]
            coeff = evaluate(control, args...; vals_dict)
            @assert coeff isa Number
            axpy!(coeff, part[1], op)
        else
            axpy!(true, part, op)
        end
    end
    return op
end


function evaluate!(op::T, generator::T, args...; kwargs...) where {T}
    if op ≡ generator
        return op
    else
        # If they're not identical, they shouldn't be of the same type. If
        # there's some weird custom type where static and timedependent
        # objects can be of the same type, they should define a custom method.
        error("typeof(op) = typeof(generator), but op ≢ generator")
    end
end


"""Substitute inside the given object.

```julia
object = substitute(object, replacements)
```

returns a modified object with the replacements defined in the given
`replacements` dictionary. Things that can be replaced include operators,
controls, and amplitudes. For example,

```julia
generator = substitute(generator::Generator, replacements)
operator = substitute(operator::Operator, replacements)
amplitude = substitute(amplitude, controls_replacements)
```

Note that `substitute` cannot be used to replace dynamic quantities, e.g.
controls, with static value. Use [`evaluate`](@ref) instead for that purpose.
"""
function substitute(object::T, replacements) where {T}
    return get(replacements, object, object)
end


function substitute(generator::Tuple, replacements::AbstractDict)
    new_generator = Any[]
    for part in generator
        if part isa Tuple
            operator, control = part
            new_operator = substitute(operator, replacements)
            new_control = substitute(control, replacements)
            push!(new_generator, (new_operator, new_control))
        else
            push!(new_generator, substitute(part, replacements))
        end
    end
    return Tuple(new_generator)
end


"""Obtain analytic parameters of the given `control`.

```julia
parameters = get_parameters(control)
```

obtains `parameters` as an `AbstractVector{Float64}` containing any tunable
analytic parameters associated with the `control`. The specific type of
`parameters` depends on how `control` is defined, but a
[`ComponentArrays.ComponentVector`](@extref) should be a common array type.

Mutating the resulting vector must directly affect the control in any
subsequent call to [`evaluate`](@ref). That is, the values in `parameters`
must alias values inside the `control`.

Note that the `control` must be an object specifically designed to have
analytic parameters. Typically, it should be implemented as a subtype of
[`ParameterizedFunction`](@ref). For a simple function `ϵ(t)` or a vector of
pulse values, which are the default types of controls discussed in the
documentation of
[`hamiltonian`](@ref QuantumPropagators.Generators.hamiltonian),
the `get_parameters` function will return an empty vector.

More generally,

```julia
parameters = get_parameters(object)
```

collects and combines all unique parameter arrays from the controls inside the
`object`. The `object` may be a
[`Generator`](@ref QuantumPropagators.Generators.Generator),
[`Trajectory`](@ref QuantumControlBase.Trajectory),
[`ControlProblem`](@ref QuantumControlBase.ControlProblem), or
any other object for which [`get_controls(object)`](@ref get_controls) is
defined. If there are multiple controls with different parameter arrays, these
are combined in a [`RecursiveArrayTools.ArrayPartition`](@extref).
This requires the `RecursiveArrayTools` package to be loaded. Again, mutating
`parameters` directly affects the underlying controls.

The `parameters` may be used as part of the `parameters`
attribute of a propagator for time-continuous dynamics, like a general ODE
solver, or in an optimization that tunes analytic control parameters, e.g.,
with a Nelder-Mead method. Examples might include the widths, peak amplitudes,
and times of a superposition of Gaussians [MachnesPRL2018](@cite), cf. the
[example of a `ParameterizedFunction`](@ref howto_parameterized),
or the amplitudes associated with spectral components in a random truncated
basis [CanevaPRA2011](@cite).

The `parameters` are not intended for optimization methods such as
[GRAPE](https://juliaquantumcontrol.github.io/GRAPE.jl/stable/) or
[Krotov](https://juliaquantumcontrol.github.io/Krotov.jl/stable/) that
fundamentally use a piecewise-constant control ansatz. In the context of such
methods, the "control parameters" are always the amplitudes of the `control` at
the mid-points of the time grid, as obtained by
[`discretize_on_midpoints`](@ref), and `get_parameters` is ignored.
"""
function get_parameters(object)
    return _get_parameters(object; via=get_controls)
end

function _get_parameters(object; via=get_controls)
    parameter_arrays = []
    seen_parameter_array = IdDict{Any,Bool}()
    seen_component = IdDict{Any,Bool}()
    for component in via(object)
        if component !== object  # E.g., get_controls(control) -> control
            if !haskey(seen_component, component)
                parameter_array = get_parameters(component)
                if isempty(parameter_array)
                    continue
                end
                if !haskey(seen_parameter_array, parameter_array)
                    push!(parameter_arrays, parameter_array)
                    seen_parameter_array[parameter_array] = true
                end
            end
            seen_component[component] = true
        end
    end
    if isempty(parameter_arrays)
        return Float64[]
    elseif length(parameter_arrays) == 1
        return parameter_arrays[1]
    else
        try
            return _combine_parameter_arrays(parameter_arrays)
        catch exception
            if exception isa MethodError
                msg = "In order for parameter arrays to be combined from multiple components, the `RecursiveArrayTools` package must be loaded"
                @error msg exception
            else
                rethrow()
            end
        end
    end
end


# Defined in extension module as
#
#   return ArrayPartition(parameter_arrays...)
#
function _combine_parameter_arrays end


"""
Abstract type for function-like objects with [`parameters`](@ref get_parameters).

A struct that is an implementation of a `ParameterizedFunction`:

* must have a `parameters` field that is an `AbstractVector` of
  floats (e.g., a [`ComponentArrays.ComponentVector`](@extref))
* must be [callable](@extref Julia Function-like-objects) with a single float
  argument `t`,
* may define getters and setters for referencing the values in `parameters`
  with convenient names.

The `parameters` field of any `ParameterizedFunction` can be accessed via
[`get_parameters`](@ref).

See [How to define a parameterized control](@ref howto_parameterized)
for an example. You may use the
[`QuantumPropagators.Interfaces.check_parameterized_function`](@ref) to check
the implementation of a `ParameterizedFunction` subtype.
"""
abstract type ParameterizedFunction <: Function end


function get_parameters(control::ParameterizedFunction)
    return getfield(control, :parameters)
end

end
