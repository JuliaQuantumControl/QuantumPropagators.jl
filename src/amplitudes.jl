module Amplitudes

export LockedAmplitude, ShapedAmplitude

import ..Controls: get_controls, evaluate, substitute, discretize_on_midpoints
using ..Controls: t_mid


#### LockedAmplitude ##########################################################


"""A time-dependent amplitude that is not a control.

```julia
ampl = LockedAmplitude(shape)
```

wraps around `shape`, which must be either a `Vector{Float64}` of values defined
on the midpoints of a time grid or a callable `shape(t)`.

```julia
ampl = LockedAmplitude(shape, tlist)
```

discretizes `shape` to the midpoints of `tlist`.
"""
struct LockedAmplitude{ST}
    shape::ST

    function LockedAmplitude(shape::ST; check = true) where {ST}
        if check && !(shape isa Vector{Float64})
            try
                shape(0.0)
            catch
                msg = "A LockedAmplitude shape must either be a Vector{Float64} or a callable"
                error(msg)
            end
        end
        return new{ST}(shape)
    end
end


function LockedAmplitude(shape, tlist)
    return LockedAmplitude(discretize_on_midpoints(shape, tlist); check = false)
end


function Base.show(io::IO, ampl::LockedAmplitude)
    print(io, "LockedAmplitude(::$(typeof(ampl.shape)))")
end


Base.Array(ampl::LockedAmplitude{Vector{Float64}}) = ampl.shape

get_controls(ampl::LockedAmplitude) = ()

function substitute(ampl::LockedAmplitude, replacements)
    return get(replacements, ampl, ampl)
end

function evaluate(ampl::LockedAmplitude{Vector{Float64}}, tlist, n::Int; _...)
    return ampl.shape[n]
end

function evaluate(ampl::LockedAmplitude{Vector{Float64}}, t::Float64; _...)
    msg = "A LockedAmplitude initialized from a vector can only be evaluated with (tlist, n)."
    error(msg)
end

function evaluate(ampl::LockedAmplitude, tlist, n::Int; _...)
    return ampl.shape(t_mid(tlist, n))
end

function evaluate(ampl::LockedAmplitude, t::Float64; _...)
    return ampl.shape(t)
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
`shape` and ``ϵ(t)`` corresponds to `control`. Each of `control` and `shape`
must be either a `Vector{Float64}` of values defined on the midpoints of a time
grid or a callable `control(t)`, respectively `shape(t)`. If both are callables,
`ampl` will also be callable. If both are vectors, they must have the same length.

```julia
ampl = ShapedAmplitude(control, tlist; shape=shape)
```

discretizes `control` and `shape` to the midpoints of `tlist`.
"""
struct ShapedAmplitude{CT,ST} <: ControlAmplitude
    control::CT
    shape::ST
end


function ShapedAmplitude(control; shape, check = true)
    if check
        if !(control isa Vector{Float64})
            try
                control(0.0)
            catch
                msg = "A ShapedAmplitude control must either be a Vector{Float64} or a callable"
                error(msg)
            end
        end
        if !(shape isa Vector{Float64})
            try
                shape(0.0)
            catch
                msg = "A ShapedAmplitude shape must either be a Vector{Float64} or a callable"
                error(msg)
            end
        end
        if (control isa Vector{Float64}) && (shape isa Vector{Float64})
            if length(control) ≠ length(shape)
                msg = "ShapedAmplitude control and shape vectors must have the same length"
                error(msg)
            end
        end
    end
    return ShapedAmplitude{typeof(control),typeof(shape)}(control, shape)
end


function Base.show(io::IO, ampl::ShapedAmplitude)
    print(io, "ShapedAmplitude(::$(typeof(ampl.control)); shape::$(typeof(ampl.shape)))")
end


function ShapedAmplitude(control, tlist::Vector{Float64}; shape)
    control = discretize_on_midpoints(control, tlist)
    shape = discretize_on_midpoints(shape, tlist)
    return ShapedAmplitude{typeof(control),typeof(shape)}(control, shape)
end


function substitute(ampl::ShapedAmplitude, replacements)
    ampl in keys(replacements) && return replacements[ampl]
    control = substitute(ampl.control, replacements)
    return ShapedAmplitude(control; shape = ampl.shape, check = true)
end


Base.Array(ampl::ShapedAmplitude{Vector{Float64},Vector{Float64}}) =
    ampl.control .* ampl.shape

(ampl::ShapedAmplitude)(t::Float64) = ampl.shape(t) * ampl.control(t)


# Vector shape: index directly; control may be vector or callable (evaluated via Controls.evaluate)
function evaluate(
    ampl::ShapedAmplitude{CT,Vector{Float64}},
    tlist::Vector{Float64},
    n::Int;
    vals_dict = IdDict()
) where {CT}
    S_t = ampl.shape[n]
    ϵ_t = evaluate(ampl.control, tlist, n; vals_dict)
    return S_t * ϵ_t
end

# Callable shape: call directly with t_mid; control may be vector or callable
function evaluate(
    ampl::ShapedAmplitude,
    tlist::Vector{Float64},
    n::Int;
    vals_dict = IdDict()
)
    S_t = ampl.shape(t_mid(tlist, n))
    ϵ_t = evaluate(ampl.control, tlist, n; vals_dict)
    return S_t * ϵ_t
end

# Callable shape and callable control: call both directly at t
function evaluate(ampl::ShapedAmplitude, t::Float64; vals_dict = IdDict())
    S_t = ampl.shape(t)
    ϵ_t = evaluate(ampl.control, t; vals_dict)
    return S_t * ϵ_t
end

function evaluate(ampl::ShapedAmplitude{Vector{Float64},Vector{Float64}}, t::Float64; _...)
    msg = "A ShapedAmplitude with vector control and shape can only be evaluated with (tlist, n)."
    error(msg)
end

function evaluate(ampl::ShapedAmplitude{Vector{Float64},ST}, t::Float64; _...) where {ST}
    msg = "A ShapedAmplitude with a vector control can only be evaluated with (tlist, n)."
    error(msg)
end

function evaluate(ampl::ShapedAmplitude{CT,Vector{Float64}}, t::Float64; _...) where {CT}
    msg = "A ShapedAmplitude with a vector shape can only be evaluated with (tlist, n)."
    error(msg)
end


function evaluate(ampl::ShapedAmplitude, args...; vals_dict = IdDict())
    S_t = evaluate(ampl.shape, args...; vals_dict)
    ϵ_t = evaluate(ampl.control, args...; vals_dict)
    return S_t * ϵ_t
end


end
