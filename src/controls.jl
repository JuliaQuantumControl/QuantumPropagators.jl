module Controls

export discretize, discretize_on_midpoints
export getcontrols
export get_tlist_midpoints
export evalcontrols, evalcontrols!, substitute_controls

using LinearAlgebra


"""Evaluate `control` at every point of `tlist`.

```julia
values = discretize(control, tlist; via_midpoints=true)
```

discretizes the given `control` to a Vector of values defined on the points of
`tlist`.

If `control` is a function, it will will first be evaluated at the midpoint of
`tlist`, see [`discretize_on_midpoints`](@ref), and then the values on the
midpoints are converted to values on `tlist`. This discretization is more
stable than directly evaluationg the control function at the values of `tlist`,
and ensures that repeated round-trips between [`discretize`](@ref) and
[`discretize_on_midpoints`](@ref) can be done safely, see the note in the
documentation of [`discretize_on_midpoints`](@ref).

The latter can still be achieved by passing `via_midpoints=false`. While such a
direct discretization is suitable e.g. for plotting, but it is unsuitable
for round-trips between [`discretize`](@ref) and
[`discretize_on_midpoints`](@ref)  (constant controls on `tlist` may result in
a zig-zag on the intervals of `tlist`).

If `control` is a vector, it will be returned un-modified if it is of the same
length as `tlist`. Otherwise, `control` must have one less value than `tlist`,
and is assumed to be defined on the midpoins of `tlist`. In that case,
[`discretize`](@ref) acts as the inverse of [`discretize_on_midpoints`](@ref).
See [`discretize_on_midpoints`](@ref) for how control values on `tlist` and
control values on the intervals of `tlist` are related.
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
        return control
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


"""Shift time grid values the interval midpoints

```julia
tlist_midpoints = get_tlist_midpoints(tlist)
```

takes a vector `tlist` of length ``n`` and returns a vector of length ``n-1``
containing the midpoint values of each interval. The intervals in `tlist` are
not required to be uniform.
"""
function get_tlist_midpoints(tlist)
    tlist_midpoints = zeros(eltype(tlist), length(tlist) - 1)
    tlist_midpoints[1] = tlist[1]
    tlist_midpoints[end] = tlist[end]
    for i = 2:length(tlist_midpoints)-1
        dt = tlist[i+1] - tlist[i]
        tlist_midpoints[i] = tlist[i] + 0.5 * dt
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
values at the beginning and end of `tlist`, in order to ensure exact bounary
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

If `control` is a vector of one less length than `tlist`, it will be returned
unchanged, under the assumption that the input is already properly discretized.

If `control` is a function, the function will be directly evaluated at the
midpoints marked as `x` in the above diagram..
"""
function discretize_on_midpoints(control::T, tlist) where {T<:Function}
    return discretize(control, get_tlist_midpoints(tlist); via_midpoints=false)
end

function discretize_on_midpoints(control::Vector, tlist)
    if length(control) == length(tlist) - 1
        return control
    elseif length(control) == length(tlist)
        vals = zeros(eltype(control), length(control) - 1)
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
controls = getcontrols(generator)
```

extracts the controls from a single dynamical generator.

For example, if `generator = hamiltonian(H0, (H1, ϵ1), (H2, ϵ2))`, extracts
`(ϵ1, ϵ2)`.
"""
getcontrols(ampl::Function) = (ampl,)
getcontrols(ampl::Vector) = (ampl,)
getcontrols(ampl::Number) = Tuple([])


"""
```julia
getcontrols(operator)
```

for a static operator (matrix) returns an empty tuple.
"""
getcontrols(operator::AbstractMatrix) = Tuple([])



"""Replace the controls in `generator` with static values.

```julia
op = evalcontrols(generator, vals_dict)
```

replaces the time-dependent controls in `generator` with the values in
`vals_dict` and returns the static operator `op`.

The `vals_dict` is a dictionary (`IdDict`) mapping controls as returned by
`getcontrols(generator)` to values.

```julia
op = evalcontrols(generator, vals_dict, tlist, n)
op = evalcontrols(generator, vals_dict, t)
```

acts similarly for generators that have explicit time dependencies (apart for
the controls). The additional parameters indicate that `generator` has explicit
time dependencies that are well-defined on the intervals of a time grid
`tlist`, and that `vals_dict` should be assumed to relate to the n'th interval
of that time grid. Respectively, an additional parameter `t` indicates the
`generator` is a time-continuous explicit dependency and that `vals_dict`
contains values defined at time `t`.

# See also:

* [`evalcontrols!`](@ref) to update `op` with new `vals_dict`.
"""
function evalcontrols(generator::Tuple, vals_dict::AbstractDict, _...)
    if isa(generator[1], Tuple)
        control = generator[1][2]
        op = vals_dict[control] * generator[1][1]
    else
        op = copy(generator[1])
    end
    for part in generator[2:end]
        if isa(part, Tuple)
            control = part[2]
            op += vals_dict[control] * part[1]
        else
            op += part
        end
    end
    return op
end

evalcontrols(num::Number, vals_dict, _...) = num

evalcontrols(operator::AbstractMatrix, _...) = operator


"""Update an existing evaluation of a `generator`.

```julia
evalcontrols!(op, generator, vals_dict, args...)
```

performs an in-place update on an `op` the was obtained from a previous call to
[`evalcontrols`](@ref) with the same `generator`, but a different `val_dict`.
"""
function evalcontrols!(op, generator::Tuple, vals_dict::AbstractDict, _...)
    if generator[1] isa Tuple
        control = generator[1][2]
        copyto!(op, generator[1][1])
        lmul!(vals_dict[control], op)
    else
        copyto!(op, generator[1])
    end
    for part in generator[2:end]
        if part isa Tuple
            control = part[2]
            axpy!(vals_dict[control], part[1], op)
        else
            axpy!(true, part, op)
        end
    end
    return op
end


"""Substitute the controls inside a `generator` with different `controls`.

```julia
new_generator = substitute_controls(generator, controls_map)
```

Creates a new generator from `generator` by replacing any control that is in
the dict `controls_map` with `controls_map[control]`. Controls that are not in
`controls_map` are kept unchanged.

The substituted controls must be time-dependent; to substitute static values
for the controls, converting the time-dependent `generator` into a static
operator, use [`evalcontrols`](@ref).

Calling `substitute_controls` on a static operator will return it unchanged.
"""
function substitute_controls(generator::Tuple, controls_map)
    new_generator = Any[]
    for part in generator
        if part isa Tuple
            operator, control = part
            new_part = (operator, get(controls_map, control, control))
            push!(new_generator, new_part)
        else
            push!(new_generator, part)
        end
    end
    return Tuple(new_generator)
end


# A static operator has no controls and remains unchanged
substitute_controls(operator::AbstractMatrix, controls_map) = operator


end
