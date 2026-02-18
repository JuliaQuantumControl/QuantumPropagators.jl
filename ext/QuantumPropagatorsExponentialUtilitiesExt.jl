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
import QuantumPropagators: init_prop, prop_step!, set_t!


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
    expv_kwargs::NamedTuple
    const timing_data::TimerOutput
end


set_t!(propagator::ExpvPropagator, t) = _pwc_set_t!(propagator, t)


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
    _...
)
```

 initializes an `ExpvPropagator`.

# Method-specific keyword arguments

* `expv_kwargs`: NamedTuple of keyword arguments forwarded to
  `ExponentialUtilities.expv`.
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
        t = float(tlist[n+1])
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
        expv_kwargs,
        timing_data,
    )
end


function prop_step!(propagator::ExpvPropagator)
    @timeit_debug propagator.timing_data "prop_step!" begin
        n = propagator.n
        tlist = getfield(propagator, :tlist)
        (0 < n < length(tlist)) || return nothing
        dt = tlist[n+1] - tlist[n]
        if propagator.backward
            dt = -dt
        end

        if propagator.inplace
            if supports_inplace(propagator.genop)
                _pwc_set_genop!(propagator, n)
                H = propagator.genop
            else
                H = _pwc_get_genop(propagator, n)
            end
            @timeit_debug propagator.timing_data "expv" begin
                Ψ = ExponentialUtilities.expv(
                    -1im * dt,
                    H,
                    propagator.state;
                    propagator.expv_kwargs...
                )
            end
            copyto!(propagator.state, Ψ)
        else
            H = _pwc_get_genop(propagator, n)
            @timeit_debug propagator.timing_data "expv" begin
                Ψ = ExponentialUtilities.expv(
                    -1im * dt,
                    H,
                    propagator.state;
                    propagator.expv_kwargs...
                )
            end
            setfield!(propagator, :state, Ψ)
        end

        _pwc_advance_time!(propagator)
        return propagator.state
    end
end

end
