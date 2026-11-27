# SPDX-FileCopyrightText: © 2021 Michael Goerz <mail@michaelgoerz.net>
#
# SPDX-License-Identifier: MIT

module Controls

export discretize, discretize_on_midpoints
export get_controls, get_parameters
export ParameterizedFunction, StaticParameter
export get_tlist_midpoints, t_mid
export evaluate, evaluate!, substitute

using LinearAlgebra
using RecursiveArrayTools: ArrayPartition

# Public, but not exported (for implementing `get_parameters` of custom types)
@static if VERSION >= v"1.11.0-DEV.469"
    eval(Meta.parse("public combine_parameters"))
end


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
function discretize(control::Function, tlist; via_midpoints = true)
    if via_midpoints
        vals_on_midpoints = discretize_on_midpoints(control, tlist)
        return discretize(vals_on_midpoints, tlist)
    else
        return Float64[control(t) for t in tlist]
    end
end

function discretize(control::Vector, tlist)
    if length(control) == length(tlist)
        return Vector{Float64}(control)
    elseif length(control) == length(tlist) - 1
        # convert `control` on intervals to values on `tlist`
        # cf. pulse_onto_tlist in Python krotov package
        vals = zeros(Float64, length(control) + 1)
        vals[1] = control[1]
        vals[end] = control[end]
        for i = 2:(length(vals)-1)
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

# See also

* [`t_mid`](@ref) – get a particular midpoint.
"""
function get_tlist_midpoints(
    tlist::AbstractVector;
    preserve_start = true,
    preserve_end = true
)
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
    for i = 2:(length(tlist_midpoints)-1)
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

# See also

* [`get_tlist_midpoints`](@ref) – get all the midpoints on which the control
  will be discretized.
* [`t_mid`](@ref) – get a particular midpoint.
* [`discretize`](@ref) – discretize directly on `tlist` instead of on the
  midpoints
"""
function discretize_on_midpoints(control::T, tlist) where {T<:Function}
    tlist_midpoints = get_tlist_midpoints(tlist)
    return discretize(control, tlist_midpoints; via_midpoints = false)
end

function discretize_on_midpoints(control::Vector, tlist)
    if length(control) == length(tlist) - 1
        return Vector{Float64}(control)
    elseif length(control) == length(tlist)
        vals = Vector{Float64}(undef, length(tlist) - 1)
        vals[1] = control[1]
        vals[end] = control[end]
        for i = 2:(length(vals)-1)
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


# Empty `vals_dict`, the default in all `evaluate` methods. Unlike a new
# `IdDict()`, it does not allocate, and lookups in it infer as "not found".
struct _EmptyValsDict <: AbstractDict{Any,Any} end
Base.length(::_EmptyValsDict) = 0
Base.isempty(::_EmptyValsDict) = true
Base.iterate(::_EmptyValsDict, state...) = nothing
Base.haskey(::_EmptyValsDict, key) = false
Base.get(::_EmptyValsDict, key, default) = default
Base.getindex(::_EmptyValsDict, key) = throw(KeyError(key))
const _NO_VALS = _EmptyValsDict()

struct _NotFound end

# Look up `key` in `vals_dict` (`_NotFound()` if missing). Checking `isempty`
# first avoids boxing an immutable `key` (e.g., a `StaticParameter`) for the
# `IdDict` lookup.
@inline _lookup(vals_dict, key) =
    isempty(vals_dict) ? _NotFound() : get(vals_dict, key, _NotFound())


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
function evaluate(object, args...; vals_dict = _NO_VALS)
    # Fallback: If `object` has to components, just look up `object` in
    # `vals_dict`
    return get(vals_dict, object, object)
end


function evaluate(operator::AbstractMatrix, args...; kwargs...)
    # objects without any explicit or explicit time dependency evaluate to
    # themselves
    return operator
end


"""
Midpoint of n'th interval of tlist.

```julia
t = t_mid(tlist, n)
```

returns the `t` that is the midpoint between points `tlist[n+1]` and
`tlist[n]`, but snapping to the beginning/end to follow the convention
explained in [`discretize_on_midpoints`](@ref) (to preserve exact boundary
conditions at the edges of the time grid.)

# See also

* [`get_tlist_midpoints`](@ref) – get all the midpoints in one go.
"""
function t_mid(tlist, n)
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


function evaluate(func::Function, tlist::Vector, n::Int64; vals_dict = _NO_VALS)
    val = _lookup(vals_dict, func)
    if val isa _NotFound
        return func(t_mid(tlist, n))
    else
        return val
    end
end


function evaluate(func::Function, t::Float64; vals_dict = _NO_VALS)
    val = _lookup(vals_dict, func)
    if val isa _NotFound
        return func(t)
    else
        return val
    end
end


function evaluate(control::Vector, tlist::Vector, n::Int64; vals_dict = _NO_VALS)
    val = _lookup(vals_dict, control)
    if !(val isa _NotFound)
        return val
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


function evaluate(control::Vector, t::Float64; vals_dict = _NO_VALS)
    val = _lookup(vals_dict, control)
    if !(val isa _NotFound)
        return val
    else
        error(
            "`evaluate(control::Vector, t::Float64)` is invalid. Use e.g. `evaluate(…, tlist, n)`."
        )
    end
end


function evaluate(generator::Tuple, args...; vals_dict = _NO_VALS)
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
function evaluate!(op, generator::Tuple, args...; vals_dict = _NO_VALS)
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
        # there's some weird custom type where static and time-dependent
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


"""Obtain the analytic parameters of an object.

```julia
parameters = get_parameters(object)
```

returns a [`RecursiveArrayTools.ArrayPartition`](@extref) of the
[parameter blocks](@ref ParameterBlocks) that `object` holds or inherits from
the controls it contains, each block once, in order of first occurrence.
Writing to `parameters` writes into the blocks, and thus changes what
[`evaluate`](@ref) returns for `object`. As long as `object` holds and
contains the same blocks, repeated calls return identical (`===`) results.

The `object` may be a control, an amplitude, a
[`Generator`](@ref QuantumPropagators.Generators.Generator),
a [`Trajectory`](@ref QuantumControl.Trajectory),
a [`ControlProblem`](@ref QuantumControl.ControlProblem), or any other object
for which [`get_controls`](@ref) is defined. A
[`ParameterizedFunction`](@ref) holds its `parameters` field as a single block.
A plain function `ϵ(t)` or a pulse holds no block. A custom type that holds
blocks in its own fields must implement its own `get_parameters` method, see
[How to implement a generator with its own parameters](@ref howto_custom_parameters).

The blocks are `parameters.x`, a tuple in the order of the index ranges that
the blocks occupy in `parameters`. Use it, e.g., to access the named fields of
a [`ComponentArrays.ComponentVector`](@extref) block.

The `parameters` are for methods that tune analytic parameters, e.g., the
widths and peak amplitudes of a superposition of Gaussians
[MachnesPRL2018](@cite), or the amplitudes of spectral components in a random
truncated basis [CanevaPRA2011](@cite). Methods that tune pulses, like
[GRAPE](@extref GRAPE :doc:`index`) and
[Krotov's method](@extref Krotov :doc:`index`), use
[`discretize_on_midpoints`](@ref) instead, and ignore `get_parameters`.

Each call goes through all controls and allocates, so call `get_parameters`
once when setting up a propagation or an optimization, and keep the result.
See [Time-Dependent Controls](@ref TimeDependentControls) for the definition of
the terms, and the how-to on
[working with control parameters](@ref howto_parameters) for examples.
"""
function get_parameters(object)
    parameters = Any[]
    seen_control = IdDict{Any,Bool}()
    for control in get_controls(object)
        if (control !== object) && !haskey(seen_control, control)
            # `control !== object`: E.g., get_controls(ϵ) -> (ϵ,)
            push!(parameters, get_parameters(control))
            seen_control[control] = true
        end
    end
    return combine_parameters(parameters...)
end


"""Combine the control `parameters` of several objects.

```julia
parameters = combine_parameters(parameters_1, parameters_2, …)
```

returns a [`RecursiveArrayTools.ArrayPartition`](@extref) of the
[parameter blocks](@ref ParameterBlocks) in the arguments, each block once, in
order of first occurrence. Each argument is the result of
[`get_parameters`](@ref) for some object, or a single block. Thus,
`combine_parameters(block).x === (block,)` for a non-empty block. An empty
block, and a block that is identical (`===`) to an earlier block, are dropped.
Two remaining blocks that may share memory, e.g., two overlapping views of the
same array, throw an `ArgumentError`. This check is not exhaustive: use
separate arrays for separate blocks, and verify with
[`check_parameterized`](@ref QuantumPropagators.Interfaces.check_parameterized).

This is a helper function for implementing [`get_parameters`](@ref) for a
custom type that holds blocks in its own fields, e.g.,

```julia
get_parameters(generator::MyGenerator) = combine_parameters(
    generator.parameters,
    (get_parameters(control) for control in get_controls(generator))...
)
```

See [How to implement a generator with its own parameters](@ref howto_custom_parameters).
"""
function combine_parameters(parameters...)
    blocks = Any[]
    for p in parameters
        # `p` is the result of `get_parameters` or a single block
        for block in ((p isa ArrayPartition) ? p.x : (p,))
            if block isa ArrayPartition
                msg = "A parameter block must not be an `ArrayPartition`. That type is reserved for the parameters of an object."
                throw(ArgumentError(msg))
            end
            if !isempty(block) && !any(b -> b === block, blocks)
                push!(blocks, block)
            end
        end
    end
    N = length(blocks)
    for i = 1:N, j = (i+1):N
        if Base.mightalias(blocks[i], blocks[j])
            msg = "Parameter blocks $i and $j (of $N) may share memory. Any two parameter blocks must either be identical (`===`) or share no memory."
            throw(ArgumentError(msg))
        end
    end
    if N == 0
        return ArrayPartition{Float64,Tuple{}}(())
    else
        return ArrayPartition(blocks...)
    end
end


"""
Abstract type for function-like objects with [`parameters`](@ref get_parameters).

A struct that is an implementation of a `ParameterizedFunction`:

* must have a `parameters` field that is a
  [parameter block](@ref ParameterBlocks) (e.g., a
  [`ComponentArrays.ComponentVector`](@extref))
* must be [callable](@extref Julia Function-like-objects) with a single float
  argument `t`,
* may define getters and setters for referencing the values in `parameters`
  with convenient names.

The `parameters` field of a `ParameterizedFunction` is the single block that
it holds: [`get_parameters`](@ref) returns an `ArrayPartition` that wraps it.
As a parameter block, it must meet four requirements:

1. It is an `AbstractVector{Float64}` with one-based indexing.
2. It can be mutated in place, and the function reads its values when it is
   evaluated, so that writing to it changes the value of `f(t)`.
3. It is not a [`RecursiveArrayTools.ArrayPartition`](@extref). That type is
   reserved for the parameters of an object with several blocks, see
   [`get_parameters`](@ref).
4. Any two blocks are either identical (`===`) or share no memory. Two
   parameterized functions share parameters by holding the same block.

See [How to define a parameterized control](@ref howto_parameterized)
for an example. You may use the
[`QuantumPropagators.Interfaces.check_parameterized_function`](@ref) to check
the implementation of a `ParameterizedFunction` subtype. See
[Time-Dependent Controls](@ref TimeDependentControls) for the rationale.
"""
abstract type ParameterizedFunction <: Function end


function get_parameters(control::ParameterizedFunction)
    return combine_parameters(getfield(control, :parameters))
end


"""A tunable parameter that is constant in time.

```julia
f = StaticParameter(parameters, i)
```

is a [`ParameterizedFunction`](@ref) that holds the
[parameter block](@ref ParameterBlocks) `parameters` and returns
`parameters[i]` for any time `t`, i.e., `f(t) == parameters[i]`. The index
`i` is optional, with a default of `1`.

This allows for tunable parameters of the drift or of a dissipator in a
[`Generator`](@ref QuantumPropagators.Generators.Generator), which does not hold
blocks itself, e.g.,

```julia
w = [ν, α]  # frequency and anharmonicity of a transmon
H = hamiltonian(
    (n̂, StaticParameter(w, 1)),
    (n̂ * (n̂ - I) / 2, StaticParameter(w, 2)),
    (a + a', ϵ)
)
```

where `w` is the only block in `get_parameters(H)` (if `ϵ` has no parameters).
"""
struct StaticParameter{P<:AbstractVector} <: ParameterizedFunction
    parameters::P
    index::Int
    function StaticParameter(parameters::P, i::Integer = 1) where {P<:AbstractVector}
        checkbounds(parameters, i)
        return new{P}(parameters, i)
    end
end

(f::StaticParameter)(t) = f.parameters[f.index]

function Base.show(io::IO, f::StaticParameter)
    print(io, "StaticParameter(")
    show(io, f.parameters)
    print(io, ", $(f.index))")
end

end
