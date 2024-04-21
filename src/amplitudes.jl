module Amplitudes

export LockedAmplitude, ShapedAmplitude

import ..Controls: get_controls, evaluate, substitute, discretize_on_midpoints
using ..Controls: t_mid


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

get_controls(ampl::LockedAmplitude) = ()

function substitute(ampl::LockedAmplitude, replacements)
    return get(replacements, ampl, ampl)
end

function evaluate(ampl::LockedPulseAmplitude, tlist, n; _...)
    return ampl.shape[n]
end

function evaluate(ampl::LockedPulseAmplitude, t; _...)
    error(
        "A LockedAmplitude initialized on a `tlist` can only be evaluated with arguments `tlist, n`."
    )
end

function evaluate(ampl::LockedContinuousAmplitude, tlist, n; _...)
    return ampl(t_mid(tlist, n))
end

function evaluate(ampl::LockedContinuousAmplitude, t; _...)
    return ampl(t)
end


#### ControlAmplitude #########################################################


# An amplitude that has `control` as the first field
abstract type ControlAmplitude end

get_controls(ampl::ControlAmplitude) = (ampl.control,)

function substitute(ampl::CT, replacements) where {CT<:ControlAmplitude}
    if ampl in keys(replacements)
        return replacements[ampl]
    else
        control = substitute(ampl.control, replacements)
        CT(control, (getfield(ampl, field) for field in fieldnames(CT)[2:end])...)
    end
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


function evaluate(ampl::ShapedAmplitude, args...; vals_dict=IdDict())
    S_t = evaluate(ampl.shape, args...; vals_dict)
    ϵ_t = evaluate(ampl.control, args...; vals_dict)
    return S_t * ϵ_t
end


end
