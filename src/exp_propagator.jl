using .Controls: get_controls
using .Interfaces: supports_inplace
using TimerOutputs: reset_timer!, @timeit_debug, TimerOutput

"""Propagator for propagation via direct exponentiation
(`method=QuantumPropagators.ExpProp`)

This is a [`PWCPropagator`](@ref).
"""
mutable struct ExpPropagator{GT,OT,ST,WT<:ExpProp.ExpPropWrk} <: PWCPropagator
    const generator::GT
    state::ST
    t::Float64  # time at which current `state` is defined
    n::Int64 # index of next interval to propagate
    const tlist::Vector{Float64}
    parameters::AbstractDict
    controls
    genop::OT
    const wrk::WT
    backward::Bool
    inplace::Bool
    convert_state::Type
    convert_operator::Type
    func::Function
    const timing_data::TimerOutput
end


set_t!(propagator::ExpPropagator, t) = _pwc_set_t!(propagator, t)


# We may want to choose the defaults for convert_state, convert_operator in
# init_prop based on the type of objects we are dealing with. See e.g.
# GradVector/GradgenOperator
_exp_prop_convert_state(state) = typeof(state)
_exp_prop_convert_operator(::Any) = Any
_exp_prop_convert_operator(::Generator) = Matrix{ComplexF64}
_exp_prop_convert_operator(::Operator) = Matrix{ComplexF64}
_exp_prop_convert_operator(::ScaledOperator) = Matrix{ComplexF64}

"""
```julia
using QuantumPropagators: ExpProp

exp_propagator = init_prop(
    state,
    generator,
    tlist;
    method=ExpProp,
    inplace=QuantumPropagators.Interfaces.supports_inplace(state),
    backward=false,
    verbose=false,
    parameters=nothing,
    func=(H_dt -> exp(-1im * H_dt))
    convert_state=_exp_prop_convert_state(state),
    convert_operator=_exp_prop_convert_operator(generator),
    _...
)
```

initializes an [`ExpPropagator`](@ref).

# Method-specific keyword arguments

* `func`: The function to evaluate. The argument `H_dt` is obtained by
  constructing an operator `H` from `generator` via the [`evaluate`](@ref)
  function and the multiplied with the time step `dt` for the current time
  interval. The propagation then simply multiplies the return value of `func`
  with the current state
* `convert_state`:  Type to which to temporarily convert states before
  multiplying the return value of `func`.
* `convert_operator`: Type to which to convert the operator `H` before
  multiplying it with `dt` and plugging the result into `func`

The `convert_state` and `convert_operator` parameters are useful for when the
`generator` and or `state` are unusual data structures for which the relevant
methods to calculate `func` are not defined. Often, it is easier to temporarily
convert them to standard complex matrices and vectors than to implement the
missing methods.
"""
function init_prop(
    state,
    generator,
    tlist,
    method::Val{:ExpProp};
    inplace = supports_inplace(state),
    backward = false,
    verbose = false,
    parameters = nothing,
    func = (H_dt -> exp(-1im * H_dt)),
    convert_state = _exp_prop_convert_state(state),
    convert_operator = _exp_prop_convert_operator(generator),
    _...
)
    if !isconcretetype(convert_state)
        throw(ArgumentError("convert_state $convert_state must be a concrete type"))
        # this is because convert_state is the `T` for `ExpPropWrk{T}`, and we
        # want `ExpPropWrk` to be a concrete type.
        # In contrast, `convert_operator` isn't used in any parametrization,
        # so it can be an abstract type.
    end
    tlist = convert(Vector{Float64}, tlist)
    controls = get_controls(generator)
    G = _pwc_get_max_genop(generator, controls, tlist)

    parameters = _pwc_process_parameters(parameters, controls, tlist)
    timing_data = TimerOutput()
    wrk = ExpProp.ExpPropWrk(convert(convert_state, state); _timing_data = timing_data)
    n = 1
    t = tlist[1]
    if backward
        n = length(tlist) - 1
        t = float(tlist[n+1])
    end
    GT = typeof(generator)
    OT = typeof(G)
    ST = typeof(state)
    WT = typeof(wrk)
    return ExpPropagator{GT,OT,ST,WT}(
        generator,
        inplace ? copy(state) : state,
        t,
        n,
        tlist,
        parameters,
        controls,
        G,
        wrk,
        backward,
        inplace,
        convert_state,
        convert_operator,
        func,
        timing_data,
    )
end


# Aliases (for backwards compatibility)
init_prop(state, generator, tlist, method::Val{:expprop}; kwargs...) =
    init_prop(state, generator, tlist, Val(:ExpProp); kwargs...)


function prop_step!(propagator::ExpPropagator)
    @timeit_debug propagator.timing_data "prop_step!" begin
        H = propagator.genop
        n = propagator.n
        tlist = getfield(propagator, :tlist)
        (0 < n < length(tlist)) || return nothing
        dt = tlist[n+1] - tlist[n]
        if propagator.backward
            dt = -dt
        end
        Ψ = convert(propagator.convert_state, propagator.state)
        if propagator.inplace
            if supports_inplace(propagator.genop)
                _pwc_set_genop!(propagator, n)
                H = convert(propagator.convert_operator, propagator.genop)
            else
                H = convert(propagator.convert_operator, _pwc_get_genop(propagator, n))
            end
            ExpProp.expprop!(Ψ, H, dt, propagator.wrk; func = propagator.func)
            if Ψ ≢ propagator.state  # `convert` of Ψ may have been a no-op
                copyto!(propagator.state, convert(typeof(propagator.state), Ψ))
            end
        else
            H = convert(propagator.convert_operator, _pwc_get_genop(propagator, n))
            Ψ = ExpProp.expprop(Ψ, H, dt, propagator.wrk; func = propagator.func)
            setfield!(propagator, :state, convert(typeof(propagator.state), Ψ))
        end
        _pwc_advance_time!(propagator)
        return propagator.state
    end
end
