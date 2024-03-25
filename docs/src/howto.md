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


## [How to define a parameterized control](@id howto_parameterized)

Parameterized controls are [function-like objects](@extref Julia Function-like-objects) with an associated vector of parameter values that must be accessible via [`QuantumPropagators.Controls.get_parameters`](@ref).

It is recommended to define a parameterized control as a subtype of [`QuantumPropagators.Controls.ParameterizedFunction`](@ref). The packages [`ComponentArrays`](https://github.com/jonniedie/ComponentArrays.jl) and [`UnPack`](https://github.com/mauro3/UnPack.jl) might be useful in the implementing of a suitable type . For example,


```jldoctest
using ComponentArrays
using UnPack: @unpack
using QuantumPropagators.Controls: ParameterizedFunction, get_parameters

struct GaussianControl <: ParameterizedFunction
    parameters::ComponentVector{Float64,Vector{Float64},Tuple{Axis{(A=1, t0=2, sigma=3)}}}
end

function GaussianControl(; A=1.0, t0=0.0, sigma=1.0)
    return GaussianControl(ComponentVector(; A, t0, sigma))
end

function (control::GaussianControl)(t)
    @unpack A, t0, sigma = control.parameters
    return A * exp(- (t - t0)^2 / (2 * sigma^2))
end

# usage

gaussian = GaussianControl(A=2.0, sigma=0.5)
gaussian.parameters.t0 = 5  # shift center from original 0.0

@show get_parameters(gaussian)
println("gaussian(4.5) = $(round(gaussian(4.5); digits=3))")

# output

get_parameters(gaussian) = (A = 2.0, t0 = 5.0, sigma = 0.5)
gaussian(4.5) = 1.213
```


We could put some extra effort into giving direct property access to all
parameters and to provide unicode-aliases for all parameters:


```jldoctest
using ComponentArrays
using QuantumPropagators.Controls: ParameterizedFunction, get_parameters

struct GaussianControl <: ParameterizedFunction
    parameters::ComponentVector{Float64,Vector{Float64},Tuple{Axis{(A=1, t0=2, sigma=3)}}}
end

function GaussianControl(; A=1.0, t0=0.0, t₀=t0, sigma=1.0, σ=sigma)
    return GaussianControl(ComponentVector(; A, t0=t₀, sigma=σ))
end

function Base.propertynames(g::GaussianControl, private::Bool=false)
    names = (:A, :t0, :t₀, :sigma, :σ)
    return private ? Tuple(union(names, fieldnames(GaussianControl))) : names
end

function Base.getproperty(g::GaussianControl, name::Symbol)
    unicode_aliases = Dict(:σ => :sigma, :t₀ => :t0)
    getproperty(get_parameters(g), get(unicode_aliases, name, name))
end

function Base.setproperty!(g::GaussianControl, name::Symbol, value)
    unicode_aliases = Dict(:σ => :sigma, :t₀ => :t0)
    setproperty!(
        get_parameters(g),
        get(unicode_aliases, name, name),
        value
    )
end

function (control::GaussianControl)(t)
    A, t₀, σ = get_parameters(control)
    return A * exp(- (t - t₀)^2 / (2 * σ^2))
end

# usage

gaussian = GaussianControl(A=2.0, σ=0.5)
gaussian.t₀ = 5  # shift center from original 0.0

@show get_parameters(gaussian)
println("gaussian(4.5) = $(round(gaussian(4.5); digits=3))")

# output

get_parameters(gaussian) = (A = 2.0, t0 = 5.0, sigma = 0.5)
gaussian(4.5) = 1.213
```

The [`QuantumPropagators.Interfaces.check_parameterized_function`](@ref) can be used to verify the implementation of a [`ParameterizedFunction`](@ref QuantumPropagators.Controls.ParameterizedFunction).
