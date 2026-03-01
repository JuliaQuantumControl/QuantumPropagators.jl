using TimerOutputs: TimerOutput


"""Propagator for Krylov `expv` propagation (`method=ExponentialUtilities`).

This is a [`PWCPropagator`](@ref). Methods that depend on
`ExponentialUtilities` are loaded via Julia package extensions.
"""
mutable struct ExponentialUtilitiesPropagator{GT,OT,ST,KST,CT} <: PWCPropagator
  const generator::GT
  state::ST
  t::Float64  # time at which current `state` is defined
  n::Int64 # index of next interval to propagate
  const tlist::Vector{Float64}
  parameters::AbstractDict
  controls
  genop::OT
  Ks::KST
  cache::CT
  backward::Bool
  inplace::Bool
  expv_kwargs::NamedTuple
  const timing_data::TimerOutput
end


set_t!(propagator::ExponentialUtilitiesPropagator, t) = _pwc_set_t!(propagator, t)
