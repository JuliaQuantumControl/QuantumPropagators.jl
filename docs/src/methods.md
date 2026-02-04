# Propagation Methods

```math
\gdef\op#1{\hat{#1}}
\gdef\ket#1{\vert{#1}\rangle}
\gdef\Liouvillian{\mathcal{L}}
\gdef\Re{\operatorname{Re}}
\gdef\Im{\operatorname{Im}}
```

As discussed in the [Overview](@ref overview_approaches), time propagation can be implemented in one of two ways:

1. By *analytically* solving the equation of motion and numerically evaluating the application time evolution operator.

   We consider this especially in the piecewise-constant case (`pwc=true` in [`propagate`](@ref)/[`init_prop`](@ref)), which is required for the traditional optimization methods [GRAPE](https://juliaquantumcontrol.github.io/GRAPE.jl/stable/) and [Krotov](https://juliaquantumcontrol.github.io/Krotov.jl/stable/). In these propagations, the time-dependent generator ``\op{H}(t)`` is [evaluated](@ref QuantumPropagators.Controls.evaluate) to a constant operator ``\op{H}`` on each interval of the time grid. The analytical solution to the Schrödinger or Liouville equation is well known, and propagation step simply has to evaluate the application of the time evolution operator ``\op{U} = \exp[-i \op{H} dt]`` to the state ``|Ψ⟩``. The following methods are built in to `QuantumPropagators`:

   * [`ExpProp`](@ref method_expprop) – constructs ``\op{U}`` explicitly and then applies it to ``|Ψ⟩``
   * [`ExponentialUtilities`](@ref method_exponentialutilities) – applies ``\op{U} |Ψ⟩`` using Krylov methods without explicitly forming ``\op{U}``
   * [`Cheby`](@ref method_cheby) — expansion of ``\op{U} |Ψ⟩`` into Chebychev polynomials, valid if ``\op{H}`` has real eigenvalues
   * [`Newton`](@ref method_newton) – expansion of ``\op{U} |Ψ⟩`` into Newton polynomials, valid if ``\op{H}`` has complex eigenvalues (non-Hermitian Hamiltonian, Liouvillian)

   The `ExpProp` method is generally not numerically efficient, but works well for small system for for debugging. The two "core" methods based on a polynomials series expansions are more suitable for bigger systems and provide both efficiency and high precision (in general, the is truncated as soon as some desired precision is reached, which is machine precision by default).

   Note that this high precision is *within the piecewise-constant approximation*. The discretization itself may introduce a non-negligible error compared to the time-continuous dynamics. There is tradeoff: A smaller `dt` decreases the discretization error, but the polynomial expansions are more effective with larger time steps.

2. By solving the equation of motion explicitly with an ODE solver.

   We support the use of any of the ODE solvers in [OrdinaryDiffEq.jl](https://docs.sciml.ai/OrdinaryDiffEq/stable/):

   * [`OrdinaryDiffEq`](@ref method_ode) – solve the equation of motion as an ODE

   The main benefit of using an ODE solver is that the generator ``\op{H}(t)`` can be treated as time-continuous, and thus avoid the time discretization error. While this is not compatible with traditional optimal control method like [GRAPE](https://juliaquantumcontrol.github.io/GRAPE.jl/stable/) and [Krotov](https://juliaquantumcontrol.github.io/Krotov.jl/stable/), it is suitable for control methods for tuning analytical control parameters [LucarelliPRA2018,MachnesPRL2018,SorensenPRA2018](@cite).

   The `method=OrdinaryDiffEq` is also available in a piecewise-constant mode by setting `pwc=true`, for comparison with `method=Cheby` and `method=Newton`.


## [`ExpProp`](@id method_expprop)

The method should be loaded with

```
using QuantumPropagators: ExpProp
```

and the passed as `method=ExpProp` to [`propagate`](@ref) or [`init_prop`](@ref):

```@docs
init_prop(state, generator, tlist, method::Val{:ExpProp}; kwargs...)
```

**Advantages**

* Simple: no knobs to turn
* "Exact" to machine precision (within the piecewise constant approximation)
* Does not require any special properties or knowledge of the dynamical generator
* Efficient for small systems

**Disadvantages**

* Bad numerical scaling with the Hilbert space dimension
* Method for `exp(-1im * H * dt)` must be defined (or `H` must be convertible to a type that can be exponentiated)

**When to use**

* Small Hilbert space dimension (<10)
* Comparing against another propagator


## [`ExponentialUtilities`](@id method_exponentialutilities)

The method requires that the [ExponentialUtilities.jl](https://docs.sciml.ai/ExponentialUtilities/stable/)
package is loaded

```
using ExponentialUtilities
```

and then passed as `method=ExponentialUtilities` (or `method=:expv`) to
[`propagate`](@ref) or [`init_prop`](@ref):

```julia
init_prop(
    state,
    generator,
    tlist,
    method::Val{:ExponentialUtilities};
    kwargs...
)
```

This method evaluates ``\exp(-i \op{H} dt) |Ψ⟩`` via a Krylov expv algorithm
without explicitly forming the matrix exponential. It is therefore often a
good fit for larger systems or matrix-free operators where direct matrix
exponentiation is too costly.

Example initialization:

```julia
using ExponentialUtilities

expv_propagator = init_prop(
    state,
    generator,
    tlist;
    method=ExponentialUtilities,  # or :expv
    inplace=QuantumPropagators.Interfaces.supports_inplace(state),
    backward=false,
    verbose=false,
    parameters=nothing,
    expv_kwargs=(; ishermitian=false),  # set for non-Hermitian generators
    convert_state=typeof(state),
    convert_operator=Matrix{ComplexF64},
)
```

**Method-specific keyword arguments**

* `expv_kwargs`: NamedTuple of keyword arguments forwarded to
  `ExponentialUtilities.expv`. Use this to set `ishermitian=false` for
  non-Hermitian generators (e.g., Liouvillians).
* `convert_state`: Type to which to temporarily convert the state before
  calling `expv`.
* `convert_operator`: Type to which to convert the operator before calling
  `expv`.

**GRAPE note**

* When using `method=:expv` with GRAPE, set `gradient_method=:taylor`.
  The default `:gradgen` path uses `GradVector`/`GradgenOperator`, which is
  not currently supported by the ExponentialUtilities Krylov backend.
* For non-Hermitian generators (e.g., Liouvillians), pass
  `prop_expv_kwargs=(; ishermitian=false)` to avoid Hermitian-only paths.

**Advantages**

* Avoids explicit construction of ``\op{U}``
* Works with matrix-free operators that support `mul!`
* Good scaling for large sparse systems

**Disadvantages**

* Requires ExponentialUtilities.jl
* Performance depends on Krylov parameters and operator structure

**When to use**

* Large, sparse, or matrix-free generators
* Piecewise-constant GRAPE-style workflows that need expv

**GRAPE note**

* When using `method=:expv` with GRAPE, set `gradient_method=:taylor`.
  The default `:gradgen` path uses `GradVector`/`GradgenOperator`, which is
  not currently supported by the ExponentialUtilities Krylov backend.
* For non-Hermitian generators (e.g., Liouvillians), pass
  `prop_expv_kwargs=(; ishermitian=false)` to avoid Hermitian-only paths.


## [`Cheby`](@id method_cheby)

The method should be loaded with

```
using QuantumPropagators: Cheby
```

and then passed as `method=Cheby` to [`propagate`](@ref) or [`init_prop`](@ref):

```@docs
init_prop(state, generator, tlist, method::Val{:Cheby}; kwargs...)
```

The time evolution operator of the piecewise-constant Schrödinger equation ``|Ψ(t)⟩ = e^{-i Ĥ dt} |Ψ(0)⟩`` is evaluated by an expansion into Chebychev polynomials [Tal-EzerJCP1984, KosloffJCP1988](@cite). This requires ``Ĥ`` to be Hermitian (have real eigenvalues) and to have a known spectral range, so that it can be normalized to the domain ``[-1, 1]`` on which the Chebychev polynomials are defined.

See [GoerzPhd2015; Chapter 3.2.1](@cite) for a detailed description of the method.

**Advantages**

* Very efficient for high precision and moderately large time steps

**Disadvantages**

* Only valid for Hermitian operators
* Must be able to estimate the spectral envelope

**When to use**

* Closed quantum systems with piecewise constant dynamics


## [`Newton`](@id method_newton)

The method should be loaded with

```
using QuantumPropagators: Newton
```

and then passed as `method=Newton` to [`propagate`](@ref) or [`init_prop`](@ref):


```@docs
init_prop(state, generator, tlist, method::Val{:Newton}; kwargs...)
```

The time evolution operator of the piecewise-constant Schrödinger equation ``|Ψ(t)⟩ = e^{-i Ĥ dt} |Ψ(0)⟩`` is evaluated by an expansion into Newton polynomials [BermanJPA1992, AshkenaziJCP1995, Tal-EzerSJSC2007](@cite). Unlike for Chebychev polynomials, this expansion does not require ``Ĥ`` to be Hermitian or to have a known spectral radius. This makes the Newton propagation applicable to open quantum systems, where ``Ĥ`` is replaced by a Liouvillian to calculate the time evolution of the density matrix.

See [GoerzPhd2015; Chapter 3.2.2](@cite) for a detailed description of the method.

**Advantages**

* Reasonably efficient for high precision and moderately large time steps
* Spectral radius does not need to be known

**Disadvantages**

* Need to choose `m_max` and `max_restarts` well for good performance.

**When to use**

* Open quantum systems with piecewise constant dynamics


## [`OrdinaryDiffEq`](@id method_ode)

The method requires that the [OrdinaryDiffEq](https://docs.sciml.ai/OrdinaryDiffEq/stable/) package is loaded

```
using OrdinaryDiffEq
```

Equivalently, the more general [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/) package can be used.

```
using DifferentialEquations
```

There is no difference between these two options: `OrdinaryDiffEq` is just a smaller dependency, but `DifferentialEquations` may be preferred if the large DifferentialEquations framework is required for the project.

In any case, the loaded package to [`propagate`](@ref) or [`init_prop`](@ref) via the `method` keyword argument:


```@docs
init_prop(state, generator, tlist, method::Val{:OrdinaryDiffEq}; kwargs...)
```

**Advantages**

* Suitable for time-continuous dynamics
* The full power of the DifferentialEquations.jl ecosystem
* Efficient for moderate precisions

**Disadvantages**

* Less efficient for piecewise-constant dynamics, and thus less suitable of PWC control methods

**When to use**

* Time-continuous dynamics
