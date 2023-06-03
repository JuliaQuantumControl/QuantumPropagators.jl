using QuantumPropagators.Controls: get_controls
using QuantumPropagators: PWCPropagator, _pwc_process_parameters, _pwc_advance_time!
using QuantumControlTestUtils.RandomObjects: random_state_vector

import QuantumPropagators: init_prop, prop_step!, set_state!, set_t!, reinit_prop!


mutable struct InvalidPropagatorNoMethods <: PWCPropagator
    generator
    state
    t::Float64  # time at which current `state` is defined
    n::Int64 # index of next interval to propagate
    tlist::Vector{Float64}
    parameters::AbstractDict
    backward::Bool
    inplace::Bool
end

function Base.show(io::IO, G::InvalidPropagatorNoMethods)
    print(io, "InvalidPropagatorNoMethods(…)")
end


function init_prop(
    state,
    generator,
    tlist,
    method::Val{:invalid_propagator_no_methods};
    inplace=true,
    backward=false,
    parameters=nothing,
    _...
)
    controls = get_controls(generator)
    parameters = _pwc_process_parameters(parameters, controls, tlist)
    n = 1
    t = tlist[1]
    if backward
        n = length(tlist) - 1
        t = float(tlist[n+1])
    end
    return InvalidPropagatorNoMethods(
        generator,
        inplace ? copy(state) : state,
        t,
        n,
        tlist,
        parameters,
        backward,
        inplace,
    )
end


"""A Propagator where the methods don't do anything"""
mutable struct InvalidPropagatorEmptyMethods <: PWCPropagator
    generator
    state
    t::Float64  # time at which current `state` is defined
    n::Int64 # index of next interval to propagate
    tlist::Vector{Float64}
    parameters
    backward::Bool
    inplace::Bool
end

function Base.show(io::IO, G::InvalidPropagatorNoMethods)
    print(io, "InvalidPropagatorNoMethods(…)")
end


function init_prop(
    state,
    generator,
    tlist,
    method::Val{:invalid_propagator_empty_methods};
    inplace=true,
    backward=false,
    parameters=nothing,
    _...
)
    n = 2 # invalid
    t = tlist[2] # invalid
    return InvalidPropagatorEmptyMethods(
        generator,
        inplace ? copy(state) : state,
        t,
        n,
        tlist,
        parameters,
        backward,
        inplace,
    )
end


function prop_step!(propagator::InvalidPropagatorEmptyMethods)
    return propagator.state
end

function set_state!(propagator::InvalidPropagatorEmptyMethods, args...)
    return nothing
end

function set_t!(propagator::InvalidPropagatorEmptyMethods, t)
    return nothing
end

function reinit_prop!(propagator::InvalidPropagatorEmptyMethods, args...; kwargs...)
    return nothing
end



"""A Propagator that spits out random wavefunctions"""
mutable struct InvalidRandomPropagator <: PWCPropagator
    generator
    state
    t::Float64  # time at which current `state` is defined
    n::Int64 # index of next interval to propagate
    tlist::Vector{Float64}
    parameters
    backward::Bool
    inplace::Bool
    seed::Int64
    rng
end


function Base.show(io::IO, G::InvalidRandomPropagator)
    print(io, "InvalidRandomPropagator(…)")
end


function init_prop(
    state,
    generator,
    tlist,
    method::Val{:invalid_random_propagator};
    inplace=true,
    backward=false,
    parameters=nothing,
    _...
)
    controls = get_controls(generator)
    parameters = _pwc_process_parameters(parameters, controls, tlist)
    n = 1
    t = tlist[1]
    if backward
        n = length(tlist) - 1
        t = float(tlist[n+1])
    end
    seed = 101134171
    rng = StableRNG(seed)
    return InvalidRandomPropagator(
        generator,
        inplace ? state : copy(state),  # the wrong way around
        t,
        n,
        tlist,
        parameters,
        backward,
        inplace,
        seed,
        rng
    )
end


function prop_step!(propagator::InvalidRandomPropagator)
    N = length(propagator.state)
    if propagator.inplace
        # mix up inpalce and not-inplace
        setfield!(propagator, :state, random_state_vector(N; propagator.rng))
    else
        copyto!(propagator.state, random_state_vector(N; propagator.rng))
    end
    return propagator.state
    # forget to call _pwc_advance_time!(propagator)
    # fail to check if the propagation is done
end


"""A Propagator forgets to return a state"""
mutable struct InvalidPropagatorNoState <: PWCPropagator
    generator
    state
    t::Float64  # time at which current `state` is defined
    n::Int64 # index of next interval to propagate
    tlist::Vector{Float64}
    parameters
    backward::Bool
    inplace::Bool
end


function Base.show(io::IO, G::InvalidPropagatorNoState)
    print(io, "InvalidPropagatorNoState(…)")
end


function init_prop(
    state,
    generator,
    tlist,
    method::Val{:invalid_propagator_no_state};
    inplace=true,
    backward=false,
    parameters=nothing,
    _...
)
    controls = get_controls(generator)
    parameters = _pwc_process_parameters(parameters, controls, tlist)
    n = 1
    t = tlist[1]
    if backward
        n = length(tlist) - 1
        t = float(tlist[n+1])
    end
    return InvalidPropagatorNoState(
        generator,
        inplace ? copy(state) : state,
        t,
        n,
        tlist,
        parameters,
        backward,
        inplace,
    )
end


function prop_step!(propagator::InvalidPropagatorNoState)
    _pwc_advance_time!(propagator)
    return nothing
end
