# Dynamical Generators

The [`propagate`](@ref QuantumPropagators.propagate) routine simulates the dynamics of a state ``|Ψ⟩`` or ``ρ̂`` under the assumption that the dynamics are described by a `generator` object. The exact equation of motion is implicit in the [`Propagator`](@ref QuantumPropagators.AbstractPropagator), but all propagators implemented in the `QuantumPropagators` package assume that the generator is the time-dependent Hamiltonian ``Ĥ(t)`` in the Schrödinger equation (``ħ=1``)

```@raw html
<a id="eq-se"></a>
```
```math
\tag{SE}
i \frac{∂}{∂t} |Ψ⟩ = Ĥ(t) |Ψ⟩\,.
```

## Operators

When evaluating the right-hand-side of Eq. [(SE)](#eq-se), the time-dependent `generator` ``Ĥ(t)`` is first evaluated into a static `operator` object ``Ĥ`` for a specific point in time via the [`QuantumPropagators.Controls.evaluate`](@ref) function. The 5-argument [`LinearAlgebra.mul!`](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#LinearAlgebra.mul!) the implements the application of the operator to the state ``|Ψ⟩``.



## Hamiltonians

The built-in [`hamiltonian`](@ref) function initializes a [`Generator`](@ref QuantumPropagators.Generators.Generator) object encapsulating a generator of the most general form

```@raw html
<a id="eq-generator"></a>
```
```math
\tag{Generator}
Ĥ(t) = Ĥ_0 + \sum_l a_l(\{ϵ_{l'}(t)\}, t) \, Ĥ_l\,.
```

The Hamiltonian consists of an (optional) drift term ``Ĥ_0`` and an arbitrary number of control terms that separate into a scalar *control amplitude* ``a_l(\{ϵ_{l'}(t)\}, t)`` and a static control operator  ``Ĥ_l``. Each *control amplitude* may consist of one or more *control function* ``ϵ_{l'}(t)``. Most commonly, ``a_l(t) ≡ ϵ_l(t)``, and thus the Hamiltonian is of the simpler form

```math
Ĥ(t) = Ĥ_0 + \sum_l ϵ_l(t) \, Ĥ_l\,.
```

The [`evaluate`](@ref QuantumPropagators.Controls.evaluate) function evaluates time-dependent [`Generator`](@ref QuantumPropagators.Generators.Generator) instances into static [`Operator`](@ref QuantumPropagators.Generators.Operator) objects.


## Control Amplitudes

The distinction between control amplitudes and control functions becomes important only in the context of optimal control, where the control functions are directly modified by optimal control, whereas the control amplitudes determine how the control functions couple to the control operators, or account for explicit time dependencies in the Hamiltonian.

Just as the [`evaluate`](@ref QuantumPropagators.Controls.evaluate) function evaluates time-dependent generators into static operators, it also evaluates control amplitudes or control functions to scalar values.


## Liouvillians

In an open quantum system, the equation of motion is assumed to take the exact same form
as in Eq. [(SE)](#eq-se),

```math
i \frac{∂}{∂t} ρ̂ = L(t) ρ̂\,,
```

where ``L`` is the Liouvillian up to a factor of ``i``.

The object representing ``L`` should be constructed with the [`liouvillian`](@ref) function, **with `convention=:TDSE`**. Just like [`hamiltonian`](@ref), this returns a [`Generator`](@ref QuantumPropagators.Generators.Generator) instance that [`evaluate`](@ref QuantumPropagators.Controls.evaluate) turns into a static [`Operator`](@ref QuantumPropagators.Generators.Operator) to be applied to a vectorized (!) state ``ρ̂``.

## Arbitrary Generators

For any Hamiltonian or Liouvillian that is more general than the form in Eq. [(Generator)](#eq-generator), a custom type would have to be implemented.  In general, the [methods defined in the `QuantumPropagators.Controls` module](@ref QuantumPropagatorsControlsAPI) (respectively [`QuantumControl.Controls`](https://juliaquantumcontrol.github.io/QuantumControl.jl/stable/api/quantum_control/#QuantumControlControlsAPI) in the broader context of optimal control) determine the relationship between generator, operators, and controls and must be implemented for that custom type.

In particular,

* [`QuantumPropagators.Controls.get_controls`](@ref) extracts the time-dependent control functions from a generator.
* [`QuantumPropagators.Controls.evaluate`](@ref) and [`evaluate!`](@ref QuantumPropagators.Controls.evaluate!) convert time-dependent generators into static operators, or control amplitudes/functions into scalar values.
* [`QuantumPropagators.Controls.substitute`](@ref) substitute control amplitudes or control functions for other control amplitudes or control functions.
