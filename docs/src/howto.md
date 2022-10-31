# Howtos

## How to extend QuantumPropagators with a new propagation method

* Define a new sub-type of [`AbstractPropagator`](@ref QuantumPropagators.AbstractPropagator) type that is unique to the propagation method, e.g. `MyNewMethodPropagator`. If appropriate, sub-type [`PiecewisePropagator`](@ref QuantumPropagators.PiecewisePropagator) or [`PWCPropagator`](@ref QuantumPropagators.PWCPropagator).

* Choose a name for the propagation method, e.g. `mynewmethod` and implement a new [`initprop`](@ref) method with the following signature

  ```
  function initprop(
      state,
      generator,
      tlist,
      method::Val{:mynewmethod};
      inplace=true,
      backward=false,
      verbose=false,
      parameters=nothing,
      # ... method-specific keyword arguments
      _...  # ignore other keyword arguments
  )
  ```

  Note the `method::Val{:mynewmethod}` as the fourth positional parameter. While the *public* interface for [`initprop`](@ref) takes `method` as a keyword argument, privately [`initprop`](@ref) dispatches for different methods as above.

* Implement the remaining methods in [The Propagator interface](@ref)
