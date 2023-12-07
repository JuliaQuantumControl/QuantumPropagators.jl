# Howtos

## How to implement a new propagation method

* Define a new sub-type of [`AbstractPropagator`](@ref QuantumPropagators.AbstractPropagator) type that is unique to the propagation method, e.g. `MyNewMethodPropagator`. If appropriate, sub-type [`PiecewisePropagator`](@ref QuantumPropagators.PiecewisePropagator) or [`PWCPropagator`](@ref QuantumPropagators.PWCPropagator).

* The high-level [`propagate`](@ref) and [`init_prop`](@ref) functions have a mandatory `method` keyword argument. That argument should receive a `Module` object for the `module` implementing a propagation method (e.g., `using QuantumPropagators: Cheby; method=Cheby`). This ensures that the module or package implementing the method is fully loaded. Internally, [`init_prop`](@ref) delegates `method=module` to a *positional* argument `method = Val(nameof(module))`

* Thus, if `MyNewMethodPropagator` is implemented in a module `MyNewMethod`, or if it wraps a registered package `MyNewMethod`, implement a new [`init_prop`](@ref) method with the following signature:

  ```
  function init_prop(
      state,
      generator,
      tlist,
      method::Val{:MyNewMethod};
      inplace=true,
      backward=false,
      verbose=false,
      parameters=nothing,
      # ... method-specific keyword arguments
      _...  # ignore other keyword arguments
  )
  ```

  Note the `method::Val{:MyNewMethod}` as the fourth positional (!) parameter. While the *public* interface for [`init_prop`](@ref) takes `method` as a keyword argument, privately [`init_prop`](@ref) dispatches for different methods as above.

* If the propagation method is not associated with a module or package, it is also possible to implement a method [`init_prop`](@ref) with a fourth positional argument with a type of, e.g., `::Val{:my_new_method}`. This would allow to call the high-level [`propagate`](@ref)/[`init_prop`](@ref) with the keyword argument `method=:my_new_method`, i.e., passing a name (`Symbol`) instead of a `Module`.

* Implement the remaining methods in [The Propagator interface](@ref).

* Test the implementation by instantiating a `propagator` and calling [`QuantumPropagators.Interfaces.check_propagator`](@ref) on it.
