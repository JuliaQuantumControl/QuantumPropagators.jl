```@meta
CurrentModule = QuantumPropagators
```

# Overview


```math
\gdef\op#1{\hat{#1}}
\gdef\ket#1{\vert{#1}\rangle}
\gdef\Liouvillian{\mathcal{L}}
\gdef\Re{\operatorname{Re}}
\gdef\Im{\operatorname{Im}}
```


The `QuantumPropagators` packages provides solvers for the dynamic equations of quantum mechanics, most importantly the Schrödinger and Liouville equations.

## Getting started

As a simple "Hello World" example, we use the [`propagate`](@ref) function to simulate a [π/2 Rabi flip in a two level system](https://gist.github.com/goerz//720669c9fffe2a8a7f4263b74d625140):

In a two-level system with ground state ``\ket{0}`` and excited state ``\ket{1}``, a constant driving field between the two levels with a pulse area of π/2 results in a population inversion, transforming the initial state ``\ket{0}`` into ``-i \ket{1}``,

```jldoctest overview
using QuantumPropagators: propagate

Ψ₀ = ComplexF64[1, 0]  #  = |0⟩
H = ComplexF64[0 1; 1 0]
tlist = collect(range(0, π/2, length=101))

Ψ = propagate(Ψ₀, H, tlist)

print("Ψ = $(round.(Ψ; digits=3))\n")

# output

Ψ = ComplexF64[0.0 + 0.0im, 0.0 - 1.0im]
```

Instead of just returning the final state, we can use `storage=true` to return an array with all the states at every point in the time grid (`tlist`). This allows us to plot the Rabi oscillation over time:

```@example overview
using Plots
gr() # hide
using QuantumPropagators: propagate # hide
Ψ₀ = ComplexF64[1, 0] # hide
H = ComplexF64[0 1; 1 0] # hide
tlist = collect(range(0, π/2, length=101)) # hide
states = propagate(Ψ₀, H, tlist; storage=true)
plot(tlist./π, abs.(states').^2; label=["ground" "excited"],
     xlabel="pulse area / π", ylabel="population", legend=:right)
```

The `storage` parameter provides a powerful way to obtain arbitrary dynamic quantities from the propagation:

* If given as `true`, return a storage array with the propagated states at each point in time instead of just the final state.
* If given a pre-allocated storage array, fill it with the propagated states at each point in time, and return the final state.
* If given in combination with `observables`, put arbitrary "observable" data derived from the propagated states in the storage array.

See [Storage of states or expectation values](@ref) for details.


## Propagation methods

When [`propagate`](@ref) is called without keyword arguments (or with `method=:auto`), an appropriate propagation method is determined automatically based on the properties of the state and generator. In the above example of a two-level system, this is `method=:expprop` which solves the Schrödinger equation by exponentiating generator to construct the time evolution operator ``\op{U} = \exp[-i \op{H} dt]`` in each time step explicitly, and applying it to the state. This is the most appropriate method for very small systems, especially a two-level system.

A specific propagation method can be forced by passing the `method` keyword argument to [`propagate`](@ref). The following methods are built-in:

* `method=:expprop`: Solve the piecewise-constant Schrödinger or Liouville equation by explicity matrix exponentiation
* `method=:cheby`: Solve the piecewise-constant Schrödinger equation (Hermitian operators)
* `method=:newton`: Solve the piecewise constant Liouville equation or non-Hermitian Schrödinger equation

The Schrödinger equation is (ħ = 1)

```math
i \frac{\partial}{\partial t} \ket{\Psi(t)} = \op{H}(t) \ket{\Psi(t)}\,.
```

For open quantum systems, we assume the Liouville equation

```math
i \frac{\partial}{\partial t} \hat{\rho}(t) = \Liouvillian(t) \hat{\rho}(t)\,,
```

which differs from most textbooks by a factor of ``i``, but has the benefit that it is structurally identical to the Schrödinger equation, so that the propagation methods do not actually need to know whether they are propagating a Hilbert space vector or a (vectorized) density matrix. See [`QuantumControl.liouvillian`](https://juliaquantumcontrol.github.io/QuantumControl.jl/stable/api/quantum_control_base/#QuantumControlBase.liouvillian) with `convention=:LvN` for how to construct an appropriate ``\Liouvillian``.



## Dynamic generators

In the above example, the "generator" `H` that is the second argument to [`propagate`](@ref) was a simple static operator. In general, we will want time-dependent Hamiltonians or Liouvillians. The standard way to initialize a time-dependent Hamiltonian is via the [`hamiltonian`](@ref) function, e.g., as  `hamiltonian(Ĥ₀, (Ĥ₁, ϵ₁), (Ĥ₂, ϵ₂))`. The `Ĥ₀`, `Ĥ₁`, and `Ĥ₂` are static operators, and `ϵ₁` and `ϵ₂` are control fields, typically functions of time `t`. For piecewise-constant propagators, `ϵ₁` nad `ϵ₂` may also be an array of amplitude values appropriate to the time grid `tlist`. The tuple-syntax for the time-dependent terms is inspired by [QuTiP](https://qutip.org/docs/latest/guide/dynamics/dynamics-time.html).

Generally, the `generator`, or the operators/controls inside the tuples can be a arbitrary objects, as long as some relevant methods are implemented for these objects, see the full section on [Dynamic Generators](@ref).

Open quantum systems are handled identically to closed quantum system, except that Hamiltonian operator are replaced by Liouvillian super-operators. For any system of non-trivial Hilbert space dimension, all (super-)operators should be sparse matrices.


## The Propagator interface

As a lower-level interface than [`propagate`](@ref), the `QuantumPropagators` package defines an interface for "propagator" objects. These are initialized via [`initprop`](@ref) as, e.g.,

```
using QuantumPropagators: initprop

propagator = initprop(Ψ₀, H, tlist)
```

The `propagator` is a propagation-method-dependent object with the interface described by [`AbstractPropagator`](@ref).

The  [`propstep!`](@ref) function can then be used to advance the `propagator`:

```@meta
DocTestSetup = quote
    using QuantumPropagators: initprop
    propagator = initprop(Ψ₀, H, tlist)
end
```

```jldoctest overview
using QuantumPropagators: propstep!

Ψ = propstep!(propagator)  # single step

while !isnothing(propstep!(propagator)); end  # go to end
Ψ = propagator.state

print("Ψ = $(round.(Ψ; digits=3)))\n")
print("t = $(round(propagator.t / π; digits=3))π\n")

# output

Ψ = ComplexF64[0.0 + 0.0im, 0.0 - 1.0im])
t = 0.5π
```

## Backward propagation

When [`propagate`](@ref) or [`initprop`](@ref) are called with `backward=true`, the propagation is initialized to run backward. The initial state is then defined at `propagator.t == tlist[end]` and each [`propstep!`](@ref) moves to the previous point in `tlist`. The equation of motion is the Schrödinger or Liouville equation with a negative $dt$. Note that the propagation uses the `generator` as it is defined: no automatic adjoint will be taken.


## Connection to DifferentialEquations.jl

The `QuantumPropagators` API is structured similarly to the [DifferentialEquations.jl](https://github.com/SciML/DifferentialEquations.jl)

* The [`propagate`](@ref) function is similar to [`DifferentialEquations.solve`](https://diffeq.sciml.ai/stable/basics/overview/#Solving-the-Problems)

* The [`initprop`](@ref) function is similar to [`DifferentialEquations.init`](https://diffeq.sciml.ai/stable/basics/integrator/#Initialization-and-Stepping)

* The [`reinitprop!`](@ref) function is similar to [`DifferentialEquations.reinit!`](https://diffeq.sciml.ai/stable/basics/integrator/#SciMLBase.reinit!)

* [The Propagator Interface](@ref) is similar to DifferentialEquations' [Integrator Interface](https://diffeq.sciml.ai/stable/basics/integrator/)

* [`propstep!`](@ref) corresponds to [`DifferentialEquations.step!`](https://diffeq.sciml.ai/stable/basics/integrator/#SciMLBase.step!)

* [`set_state!`](@ref) corresponds to [`DifferentialEquations.set_u!`](https://diffeq.sciml.ai/stable/basics/integrator/#SciMLBase.set_u!)

* [`set_t!`](@ref) corresponds to [`DifferentialEquations.set_t!`](https://diffeq.sciml.ai/stable/basics/integrator/#SciMLBase.set_t!)


Note that the equation of motion for `QuantumPropagators` is implicit in the propagation `method` (usually the Schrödinger/Liouville equation), so the initialization of a Propagator via the initial state and the "generator" is more specialized than DifferentialEquations' [Problem Interface](https://diffeq.sciml.ai/stable/basics/problem/#Problem-Interface).
