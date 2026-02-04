module QuantumPropagatorsExponentialUtilitiesExt

using LinearAlgebra
using TimerOutputs: @timeit_debug, TimerOutput
using ExponentialUtilities

using QuantumPropagators:
    QuantumPropagators,
    PWCPropagator,
    _pwc_get_max_genop,
    _pwc_get_genop,
    _pwc_set_genop!,
    _pwc_set_t!,
    _pwc_advance_time!,
    _pwc_process_parameters
using QuantumPropagators.Controls: get_controls
using QuantumPropagators.Interfaces: supports_inplace
using QuantumPropagators.Generators: ScaledOperator
import QuantumPropagators: init_prop, prop_step!, set_t!

struct _SizedOp{A}
    A::A
end

Base.size(op::_SizedOp) = size(op.A)
Base.size(op::_SizedOp, dim::Integer) = size(op.A)[dim]
Base.eltype(op::_SizedOp) = eltype(op.A)
LinearAlgebra.ishermitian(op::_SizedOp) = LinearAlgebra.ishermitian(op.A)
LinearAlgebra.mul!(y, op::_SizedOp, x) = mul!(y, op.A, x)

_ensure_size_dim(A) =
    hasmethod(size, Tuple{typeof(A), Int}) ? A :
    (hasmethod(size, Tuple{typeof(A)}) ? _SizedOp(A) : A)


"""Propagator for Krylov expv propagation via ExponentialUtilities (`method=ExponentialUtilities`).

This is a [`PWCPropagator`](@ref).
"""
mutable struct ExpvPropagator{GT,OT,ST} <: PWCPropagator
    const generator::GT
    state::ST
    t::Float64  # time at which current `state` is defined
    n::Int64 # index of next interval to propagate
    const tlist::Vector{Float64}
    parameters::AbstractDict
    controls
    genop::OT
    backward::Bool
    inplace::Bool
    convert_state::Type
    convert_operator::Type
    expv_kwargs::NamedTuple
    const timing_data::TimerOutput
end


set_t!(propagator::ExpvPropagator, t) = _pwc_set_t!(propagator, t)


_expv_convert_state(state) = typeof(state)
_expv_convert_operator(::Any) = Any
_expv_convert_operator(::QuantumPropagators.Generator) = Matrix{ComplexF64}
_expv_convert_operator(::QuantumPropagators.Operator) = Matrix{ComplexF64}
_expv_convert_operator(::QuantumPropagators.ScaledOperator) = Matrix{ComplexF64}


"""
```julia
using ExponentialUtilities

expv_propagator = init_prop(
    state,
    generator,
    tlist;
    method=ExponentialUtilities,
    inplace=QuantumPropagators.Interfaces.supports_inplace(state),
    backward=false,
    verbose=false,
    parameters=nothing,
    expv_kwargs=(;),
    convert_state=_expv_convert_state(state),
    convert_operator=_expv_convert_operator(generator),
    _...
)
```

initializes an [`ExpvPropagator`](@ref).

# Method-specific keyword arguments

* `expv_kwargs`: NamedTuple of keyword arguments forwarded to
  `ExponentialUtilities.expv`.
* `convert_state`: Type to which to temporarily convert the state before
  calling `expv`.
* `convert_operator`: Type to which to convert the operator before calling
  `expv`.
"""
function init_prop(
    state,
    generator,
    tlist,
    method::Val{:ExponentialUtilities};
    inplace = supports_inplace(state),
    backward = false,
    verbose = false,
    parameters = nothing,
    expv_kwargs = (;),
    convert_state = _expv_convert_state(state),
    convert_operator = _expv_convert_operator(generator),
    _...
)
    tlist = convert(Vector{Float64}, tlist)
    controls = get_controls(generator)
    G = _pwc_get_max_genop(generator, controls, tlist)

    parameters = _pwc_process_parameters(parameters, controls, tlist)
    timing_data = TimerOutput()

    n = 1
    t = tlist[1]
    if backward
        n = length(tlist) - 1
        t = float(tlist[n + 1])
    end
    GT = typeof(generator)
    OT = typeof(G)
    ST = typeof(state)
    return ExpvPropagator{GT,OT,ST}(
        generator,
        inplace ? copy(state) : state,
        t,
        n,
        tlist,
        parameters,
        controls,
        G,
        backward,
        inplace,
        convert_state,
        convert_operator,
        expv_kwargs,
        timing_data,
    )
end


# Aliases
init_prop(state, generator, tlist, method::Val{:expv}; kwargs...) =
    init_prop(state, generator, tlist, Val(:ExponentialUtilities); kwargs...)


function prop_step!(propagator::ExpvPropagator)
    @timeit_debug propagator.timing_data "prop_step!" begin
        if nameof(typeof(propagator.state)) == :GradVector &&
           nameof(parentmodule(typeof(propagator.state))) == :QuantumGradientGenerators
            throw(ArgumentError(
                "ExponentialUtilities propagation does not support GRAPE `gradient_method=:gradgen`. " *
                "Use `gradient_method=:taylor` instead."
            ))
        end
        H = propagator.genop
        n = propagator.n
        tlist = getfield(propagator, :tlist)
        (0 < n < length(tlist)) || return nothing
        dt = tlist[n + 1] - tlist[n]
        if propagator.backward
            dt = -dt
        end
        dt_expv = complex(dt)

        Ψ = convert(propagator.convert_state, propagator.state)
        if propagator.inplace
            if supports_inplace(propagator.genop)
                _pwc_set_genop!(propagator, n)
                H = convert(propagator.convert_operator, propagator.genop)
            else
                H = convert(propagator.convert_operator, _pwc_get_genop(propagator, n))
            end
            H = _ensure_size_dim(H)
            H = ScaledOperator(-1im, H)
            @timeit_debug propagator.timing_data "expv" begin
                Ψ = ExponentialUtilities.expv(dt_expv, H, Ψ; propagator.expv_kwargs...)
            end
            copyto!(propagator.state, convert(typeof(propagator.state), Ψ))
        else
            H = convert(propagator.convert_operator, _pwc_get_genop(propagator, n))
            H = _ensure_size_dim(H)
            H = ScaledOperator(-1im, H)
            @timeit_debug propagator.timing_data "expv" begin
                Ψ = ExponentialUtilities.expv(dt_expv, H, Ψ; propagator.expv_kwargs...)
            end
            setfield!(propagator, :state, convert(typeof(propagator.state), Ψ))
        end

        _pwc_advance_time!(propagator)
        return propagator.state
    end
end

end
