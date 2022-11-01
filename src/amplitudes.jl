module Amplitudes

export LockedAmplitude, ShapedAmplitude

import ..Controls: getcontrols, evalcontrols, substitute_controls, discretize_on_midpoints


#### LockedAmplitude ##########################################################


"""A time-dependent amplitude that is not a control.

```julia
ampl = LockedAmplitude(shape)
```

wraps around `shape`, which must be either a vector of values defined on the
midpoints of a time grid or a callable `shape(t)`.

```julia
ampl = LockedAmplitude(shape, tlist)
```

discretizes `shape` to the midpoints of `tlist`.
"""
abstract type LockedAmplitude end


function LockedAmplitude(shape)
    if shape isa Vector{Float64}
        return LockedPulseAmplitude(shape)
    else
        return LockedContinuousAmplitude(shape)
    end
end


function LockedAmplitude(shape, tlist)
    return LockedPulseAmplitude(discretize_on_midpoints(shape, tlist))
end


function Base.show(io::IO, ampl::LockedAmplitude)
    print(io, "LockedAmplitude(::$(typeof(ampl.shape)))")
end


struct LockedPulseAmplitude <: LockedAmplitude
    shape::Vector{Float64}
end


Base.Array(ampl::LockedPulseAmplitude) = ampl.shape


struct LockedContinuousAmplitude <: LockedAmplitude

    shape

    function LockedContinuousAmplitude(shape)
        try
            S_t = shape(0.0)
        catch
            error("A LockedAmplitude shape must either be a vector of values or a callable")
        end
        return new(shape)
    end

end

(ampl::LockedContinuousAmplitude)(t::Float64) = ampl.shape(t)

getcontrols(ampl::LockedAmplitude) = ()

substitute_controls(ampl::LockedAmplitude, controls_map) = ampl

function evalcontrols(ampl::LockedPulseAmplitude, vals_dict, tlist, n)
    return ampl.shape[n]
end

function evalcontrols(ampl::LockedContinuousAmplitude, vals_dict, tlist, n)
    # It's technically possible to determine t from (tlist, n), but maybe we
    # shouldn't
    error("LockedAmplitude must be initialized with tlist")
end


#### ControlAmplitude #########################################################


# An amplitude that has `control` as the first field
abstract type ControlAmplitude end

getcontrols(ampl::ControlAmplitude) = (ampl.control,)

function substitute_controls(ampl::CT, controls_map) where {CT<:ControlAmplitude}
    control = get(controls_map, ampl.control, ampl.control)
    CT(control, (getfield(ampl, field) for field in fieldnames(CT)[2:end])...)
end


#### ShapedAmplitude ##########################################################


"""Product of a fixed shape and a control.

```julia
ampl = ShapedAmplitude(control; shape=shape)
```

produces an amplitude ``a(t) = S(t) ϵ(t)``, where ``S(t)`` corresponds to
`shape` and ``ϵ(t)`` corresponds to `control`. Both `control` and `shape`
should be either a vector of values defined on the midpoints of a time grid or
a callable `control(t)`, respectively `shape(t)`. In the latter case, `ampl`
will also be callable.

```julia
ampl = ShapedAmplitude(control, tlist; shape=shape)
```

discretizes `control` and `shape` to the midpoints of `tlist`.
"""
abstract type ShapedAmplitude <: ControlAmplitude end

function ShapedAmplitude(control; shape)
    if (control isa Vector{Float64}) && (shape isa Vector{Float64})
        return ShapedPulseAmplitude(control, shape)
    else
        try
            ϵ_t = control(0.0)
        catch
            error(
                "A ShapedAmplitude control must either be a vector of values or a callable"
            )
        end
        try
            S_t = shape(0.0)
        catch
            error("A ShapedAmplitude shape must either be a vector of values or a callable")
        end
        return ShapedContinuousAmplitude(control, shape)
    end
end

function Base.show(io::IO, ampl::ShapedAmplitude)
    print(io, "ShapedAmplitude(::$(typeof(ampl.control)); shape::$(typeof(ampl.shape)))")
end

function ShapedAmplitude(control, tlist; shape)
    control = discretize_on_midpoints(control, tlist)
    shape = discretize_on_midpoints(shape, tlist)
    return ShapedPulseAmplitude(control, shape)
end

struct ShapedPulseAmplitude <: ShapedAmplitude
    control::Vector{Float64}
    shape::Vector{Float64}
end

Base.Array(ampl::ShapedPulseAmplitude) = ampl.control .* ampl.shape


struct ShapedContinuousAmplitude <: ShapedAmplitude
    control
    shape
end

(ampl::ShapedContinuousAmplitude)(t::Float64) = ampl.shape(t) * ampl.control(t)


"""
```julia
aₙ = evalcontrols(ampl, vals_dict)
```

evaluates a general control amplitude by replacing each control with the
values in `vals_dict`. Returns a number.

Note that for "trivial" amplitudes (where the amplitude is identical to the
control), this simply looks up the control in `vals_dict`.

For amplitudes with explicit time dependencies (outside of the controls),
additional parameters `tlist, n` or `t` should be given to indicate the time at
which the amplitude is to be evaluated.
"""
evalcontrols(ampl::Function, vals_dict, _...) = vals_dict[ampl]
evalcontrols(ampl::Vector, vals_dict, _...) = vals_dict[ampl]


function evalcontrols(ampl::ShapedPulseAmplitude, vals_dict, tlist, n)
    return ampl.shape[n] * vals_dict[ampl.control]
end

function evalcontrols(ampl::ShapedContinuousAmplitude, vals_dict, tlist, n)
    # It's technically possible to determind t from (tlist, n), but maybe we
    # shouldn't
    error("ShapedAmplitude must be initialized with tlist")
end


end
