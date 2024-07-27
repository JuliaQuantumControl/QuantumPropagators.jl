# Overview


```math
\gdef\op#1{\hat{#1}}
\gdef\Liouvillian{\mathcal{L}}
\gdef\Re{\operatorname{Re}}
\gdef\Im{\operatorname{Im}}
```


The `QuantumPropagators` packages provides solvers for the dynamic equations of quantum mechanics, most importantly the Schrödinger and Liouville equations. We refer to the numerical evaluation of a quantum state ``|Ψ(t + dt)⟩``  from a state ``|Ψ(t)⟩`` under a given equation of motion as "time propagation".

## [Getting started](@id getting_started)

As a simple "Hello World" example, we use the [`propagate`](@ref) function to simulate a [π/2 Rabi flip in a two level system](https://gist.github.com/goerz//720669c9fffe2a8a7f4263b74d625140):

In a two-level system with ground state ``|0⟩`` and excited state ``|1⟩``, a constant driving field between the two levels with a pulse area of π/2 results in a population inversion, transforming the initial state ``|0⟩`` into ``-i |1⟩``,

```jldoctest overview
using QuantumPropagators: propagate, ExpProp

Ψ₀ = ComplexF64[1, 0]  #  = |0⟩
H = ComplexF64[0 1; 1 0]
tlist = collect(range(0, π/2, length=101))

Ψ = propagate(Ψ₀, H, tlist; method=ExpProp)

print("Ψ = $(round.(Ψ; digits=3))\n")

# output

Ψ = ComplexF64[0.0 + 0.0im, 0.0 - 1.0im]
```

We've used a simple [exponential propagator (`ExpProp`)](@ref method_expprop) here, which directly calculates `U = exp(-1im * H * dt)` in every time step and applies it to the current state `Ψ`.

Instead of just returning the final state, we can use `storage=true` to return an array with all the states at every point in the time grid (`tlist`). This allows us to plot the Rabi oscillation over time:

```@example overview
using Plots
gr() # hide
using QuantumPropagators: propagate, ExpProp # hide
Ψ₀ = ComplexF64[1, 0] # hide
H = ComplexF64[0 1; 1 0] # hide
tlist = collect(range(0, π/2, length=101)) # hide
states = propagate(Ψ₀, H, tlist; method=ExpProp, storage=true)
plot(tlist./π, abs.(states').^2; label=["ground" "excited"],
     xlabel="pulse area / π", ylabel="population", legend=:right)
```

The `storage` parameter provides a powerful way to obtain arbitrary dynamic quantities from the propagation:

* If given as `true`, return a storage array with the propagated states at each point in time instead of just the final state.
* If given a pre-allocated storage array, fill it with the propagated states at each point in time, and return the final state.
* If given in combination with `observables`, put arbitrary "observable" data derived from the propagated states in the storage array.

See the discussion of [Expectation Values](@ref) for details.

## [Approaches to time propagation](@id overview_approaches)

We are primarily interested in the time propagation of a quantum state under
the Schrödinger equation (ħ = 1)

```math
i \frac{\partial}{\partial t} |\Psi(t)⟩ = \op{H}(t) |\Psi(t)⟩\,,
```

which describes the dynamics of a closed quantum system. For open quantum systems, the equivalent equation is the Liouville equation, which we write as

```math
i \frac{\partial}{\partial t} \hat{\rho}(t) = \Liouvillian(t) \hat{\rho}(t)\,.
```

This form differs from most textbooks by a factor of ``i``, but has the benefit that it is structurally identical to the Schrödinger equation, so that the propagation methods do not actually need to know whether they are propagating a Hilbert space vector or a (vectorized) density matrix. See [`liouvillian`](@ref) with `convention=:LvN` for how to construct an appropriate ``\Liouvillian``.

There are two fundamental approaches to solving the Schrödinger equation (or any equation of motion):

1. We can analytically solve the Schrödinger equation and then numerically evaluate the solution. Mathematically, this is the application of the time evolution operator ``\op{H}`` as ``|Ψ(t+dt)⟩ = \op{U}(t) |Ψ(t)⟩``. For a piecewise-constant ``\op{H}(t)``where there is a time-independent ``\op{H}`` in the interval ``[t, t+dt]``, the time evolution operator is well-known to be

   ```math
   \op{U} = \exp[-i \op{H} dt]\,.
   ```

   The propagated state ``\op{U} |Ψ(t)⟩`` would then be obtained, e.g., by expanding the exponentiation of the operator ``\op{H}`` into a polynomial series. This can be done to arbitrary precision by truncating the series at an appropriate point.

2.  We can use a general ODE solver, e.g., using some kind of Runge-Kutta scheme. These work by following the derivative between ``t`` and ``t+dt`` with some adaptive internal step size to achieve a given precision.

The first case of a piecewise-constant time evolution operator is particularly relevant to quantum control, since the two most venerable methods of quantum control ([GRAPE](https://juliaquantumcontrol.github.io/GRAPE.jl/stable/) and [Krotov's method](https://juliaquantumcontrol.github.io/Krotov.jl/stable/)) are inherently piecewise-constant. Hence, the `QuantumPropagators` package implements two efficient polynomial propagators that have a long history in quantum control, using Chebychev polynomials for closed quantum systems [Tal-EzerJCP1984, KosloffJCP1988](@cite) and Newton polynomials for open quantum systems [BermanJPA1992, KosloffARPC94,AshkenaziJCP1995](@cite).

For propagation via an ODE Solver, `QuantumPropagators` delegates to the [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/), respectively the [OrdinaryDiffEq.jl](https://docs.sciml.ai/OrdinaryDiffEq/stable/) package.


## [Propagation methods](@id overview_propagation)

The [`propagate`](@ref) function has a mandatory `method` keyword argument.
It should be passed a module that implements the method. For the built-in methods, this would be one of the following submodule of `QuantumPropagators`:

```julia
using QuantumPropagators: ExpProp, Cheby, Newton
```

The equation of motion is implicit in the propagation method. The above methods target the Schrödinger or Liouville equations for a piecewise-constant Hamiltonian that is evaluated on the midpoints of the propagation time grid.

The two core method are:

* `method=Cheby`: Evaluate the application of the unitary time-evolution operator via an expansion into Chebychev polynomials.
* `method=Newton`: Evaluate the application of the non-unitary time-evolution operator via an expansion into Newton polynomials, for a Liouvillian or a non-Hermitian Hamiltonian.

Furthermore, as a fallback for very small system or for debugging,

* `method=ExpProp`: Explicitly construct the time evolution operator by matrix exponentiation and apply it to the state.

If the [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) or (equivalently) th [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/) package is loaded, `QuantumPropagators` can delegate to it:

```julia
using OrdinaryDiffEq
```

allows to pass

* `method=OrdinaryDiffEq`: Use [OrdinaryDiffEq.jl](https://github.com/SciML/OrdinaryDiffEq.jl) as a backend with any of the [algorithms available for ODEs in DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/solvers/ode_solve/), `alg=Tsit5()` by default.

Unlike any of the built-in methods, `OrdinaryDiffEq` is able to propagate for time-continuous generators. This is the default for that propagator (`pwc=false`). By setting `pwc=true` or `piecewise=true`) the ODE solvers can also be used for piecewise-constant Hamiltonians or Liouvillians, providing an alternative to the built-in `method=Cheby` and `method=Newton`.

See the more extended discussion of [Propagation Methods](@ref) for more details.


## [Dynamical generators](@id overview_dynamical_generators)

In the [initial example](@ref getting_started), the "generator" `H` that is the second argument to [`propagate`](@ref) was a simple static operator. In general, we will want time-dependent Hamiltonians or Liouvillians. The standard way to initialize a time-dependent Hamiltonian is via the [`hamiltonian`](@ref) function, e.g., as  `hamiltonian(Ĥ₀, (Ĥ₁, ϵ₁), (Ĥ₂, ϵ₂))`. The `Ĥ₀`, `Ĥ₁`, and `Ĥ₂` are static operators, and `ϵ₁` and `ϵ₂` are control fields, typically functions ([or function-like objects](https://docs.julialang.org/en/v1/manual/methods/#Function-like-objects)) of time `t`. For piecewise-constant propagators, `ϵ₁` and `ϵ₂` may also be an array of amplitude values appropriate to the time grid `tlist`. The tuple-syntax for the time-dependent terms is inspired by [QuTiP](@extref qutip :label:`time`).

Generally, the `generator`, or the operators/controls inside the tuples can be a arbitrary objects, as long as some relevant methods are implemented for these objects, see the full section on [Dynamical Generators](@ref).

Open quantum systems are handled identically to closed quantum system, except that Hamiltonian operator are replaced by Liouvillian super-operators. For any system of non-trivial Hilbert space dimension, all (super-)operators should be sparse matrices.


## The Propagator interface

As a lower-level interface than [`propagate`](@ref), the `QuantumPropagators` package defines an interface for "propagator" objects. These are initialized via [`init_prop`](@ref) as, e.g.,

```
using QuantumPropagators: init_prop

propagator = init_prop(Ψ₀, H, tlist; method)
```

with a mandatory `method` keyword argument.

The `propagator` is a propagation-method-dependent object with the interface described by [`AbstractPropagator`](@ref) and [`QuantumPropagators.Interfaces.check_propagator`](@ref).

The  [`prop_step!`](@ref) function can then be used to advance the `propagator`:

```@meta
DocTestSetup = quote
    using QuantumPropagators: init_prop, ExpProp
    propagator = init_prop(Ψ₀, H, tlist; method=ExpProp)
end
```

```jldoctest overview
using QuantumPropagators: prop_step!

Ψ = prop_step!(propagator)  # single step

while !isnothing(prop_step!(propagator)); end  # go to end
Ψ = propagator.state

print("Ψ = $(round.(Ψ; digits=3)))\n")
print("t = $(round(propagator.t / π; digits=3))π\n")

# output

Ψ = ComplexF64[0.0 + 0.0im, 0.0 - 1.0im])
t = 0.5π
```

## In-place propagation

Most propagators support both an in-place and a not-in-place mode. These modes can be switched via the `inplace` parameter to [`propagate`](@ref)/[`init_prop`](@ref), which defaults to the value of [`QuantumPropagators.Interfaces.supports_inplace`](@ref) for the given initial `state`. When using in-place operations, propagators should minimize the allocation of memory and modify states using in-place linear algebra operations ([BLAS](@extref Julia `BLAS-functions`)). Otherwise, with `inplace=false`, mutating states or operators is avoided. See the in-place and not-in-place operations described in [`QuantumPropagators.Interfaces.check_state`](@ref) and [`QuantumPropagators.Interfaces.check_operator`](@ref).

In-place operations can be dramatically more efficient for large Hilbert space dimensions. On the other hand, not-in-place operations can be more efficient for small Hilbert spaces, in particular when a [static vector](@extref StaticArrays `SVector`) can be used to represent the state. Moreover, frameworks for automatic differentiation such as [Zygote](https://fluxml.ai/Zygote.jl/stable/) do not support in-place operations.

When using custom structs for states, operators, or generators, the struct itself not not need to be mutable (according to [`Base.ismutable`](@extref Julia)) in order to support `inplace=true`. It only must support the in-place operations defined in the [formal interface](@ref QuantumPropagatorsInterfacesAPI) and indicate that support by defining [`QuantumPropagators.Interfaces.supports_inplace`](@ref).
Typically, in-place operations on immutable custom structs involve mutating the mutable properties of that struct.


## Backward propagation

When [`propagate`](@ref) or [`init_prop`](@ref) are called with `backward=true`, the propagation is initialized to run backward. The initial state is then defined at `propagator.t == tlist[end]` and each [`prop_step!`](@ref) moves to the previous point in `tlist`. The equation of motion is the Schrödinger or Liouville equation with a negative $dt$. For a Hermitian `generator`, doing a forward propagation followed by a backward propagation will recover the initial state. For a non-Hermitian `generator`, this no longer holds. Note that in optimal control methods such as GRAPE or Krotov's method, obtaining gradients involves a "backward propagation with the adjoint generator" (when the generator is non-Hermitian and adjoint/non-adjoint makes a difference). The [`propagate`](@ref) routine with `backward=true` will not automatically take this adjoint of the `generator`; instead, the adjoint generator must be passed explicitly.

## Parameterized controls

Controls may depend on a list of tunable parameters. Such controls must be especially defined and should be subtypes of [`QuantumPropagators.Controls.ParameterizedFunction`](@ref), or more generally ["functors"](https://docs.julialang.org/en/v1/manual/methods/#Function-like-objects) of a single float `t` for which [`QuantumPropagators.Controls.get_parameters`](@ref)
is implemented.

Independently, all [propagators](@ref QuantumPropagators.AbstractPropagator) have a field `parameters` that is a dict of controls to *propagation parameters* for that control. For continuous-time propagators, such as a propagator initialized with `method=OrdinaryDiffEq` and `pwc=false`, these propagation parameters are exactly the analytic parameters of the `control` as obtained by [`QuantumPropagators.Controls.get_parameters`](@ref).

For piecewise-constant propagators (all the default built-in propagators), the propagation parameters are always the values of the control evaluated on the mid points of the time grid, see [`QuantumPropagators.Controls.discretize_on_midpoints`](@ref). Specifically, the analytic parameters from `get_parameters(control)` are not used as propagation parameters for piecewise propagators.

In any case, mutating `propagator.parameters` affects the propagation in subsequent calls to [`prop_step!`](@ref).


## Connection to DifferentialEquations.jl

The `QuantumPropagators` API is structured similarly to the [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl)

* The [`propagate`](@ref) function is similar to [`DifferentialEquations.solve`](https://diffeq.sciml.ai/stable/basics/overview/#Solving-the-Problems)

* The [`init_prop`](@ref) function is similar to [`DifferentialEquations.init`](https://diffeq.sciml.ai/stable/basics/integrator/#Initialization-and-Stepping)

* The [`reinit_prop!`](@ref) function is similar to [`DifferentialEquations.reinit!`](https://diffeq.sciml.ai/stable/basics/integrator/#SciMLBase.reinit!)

* [The Propagator interface](@ref) is similar to DifferentialEquations' [Integrator Interface](https://diffeq.sciml.ai/stable/basics/integrator/)

* [`prop_step!`](@ref) corresponds to [`DifferentialEquations.step!`](https://diffeq.sciml.ai/stable/basics/integrator/#SciMLBase.step!)

* [`set_state!`](@ref) corresponds to [`DifferentialEquations.set_u!`](https://diffeq.sciml.ai/stable/basics/integrator/#SciMLBase.set_u!)

* [`set_t!`](@ref) corresponds to [`DifferentialEquations.set_t!`](https://diffeq.sciml.ai/stable/basics/integrator/#SciMLBase.set_t!)


Note that the equation of motion for `QuantumPropagators` is implicit in the propagation `method` (usually the Schrödinger/Liouville equation), so the initialization of a Propagator via the initial state and the "generator" is more specialized than DifferentialEquations' [Problem Interface](https://diffeq.sciml.ai/stable/basics/problem/#Problem-Interface).

The propagator returned by `using DifferentialEquations; init_prop(Ψ₀, H, tlist; method=DifferentialEquations)` is a thin wrapper around DifferentialEquations' integrator. That propagator uses an un-exported function [`QuantumPropagators.ode_function`](@ref) to wrap around the evaluation of a time-dependent generator. The [`ode_function`](@ref QuantumPropagators.ode_function) wrapper could also be used directly to enable working with the data structures defined in `QuantumPropagators` in the context of the `DifferentialEquations` package.
