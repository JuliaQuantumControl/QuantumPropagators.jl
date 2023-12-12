# Howtos

```math
\gdef\op#1{\hat{#1}}
\gdef\Liouvillian{\mathcal{L}}
\gdef\Re{\operatorname{Re}}
\gdef\Im{\operatorname{Im}}
```

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


## How to specify the spectral range for a Chebychev propagation

A propagation with [`method=Cheby`](@ref method_cheby) requires that the dynamic generator ``\op{H}(t)`` be normalized to a spectral range of [-1, 1]. That is, the method needs a (pessimistic) estimate of the "spectral envelope": the minimum and maximum eigenvalue of ``\op{H}(t)`` for any point `t` on the interval of the propagation time grid `tlist`.

By default, the Chebychev propagator uses heuristics to estimate  this spectral envelope. If the spectral envelope is known (either analytically of via a separate numerical exploration of the eigenvalues over the full range of possible controls), the minimum and maximum eigenvalues of ``\op{H}(t)`` can be passed as keyword arguments `E_min` and `E_max` to [`propagate`](@ref) or [`init_prop`](@ref). Since the Chebychev method is only defined for Hermitian generators, `E_min` and `E_max` must be real values. Both values must be given.

Manually specifying `E_min` and `E_max` works with the default `specrange_method=:auto` as well as with the explicit `specrange_method=:manual`. When calculating the Chebychev coefficients, the given values may still be enlarged by the default `specrange_buffer` keyword argument in [`init_prop`](@ref). If `E_min` and `E_max` should be used *exactly*, pass `specrange_buffer=0`.
