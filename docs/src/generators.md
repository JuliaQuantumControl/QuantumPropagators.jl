# Dynamic Generators

The [`propagate`](@ref QuantumPropagators.propagate) routine simulates the dynamics of a state ``|Ψ⟩`` or ``ρ̂`` under the assumption that that time derivative of the state is described by a `generator`. The exact equation of motion is implicit in the [`Propagator`](@ref QuantumPropagators.AbstractPropagator), but all propagators implemented in the `QuantumPropagators` package assume that the generator is the time-dependent Hamiltonian ``Ĥ(t)`` in the Schrödinger equation (``ħ=1``)

```math
i \frac{∂}{∂t} |Ψ⟩ = Ĥ(t) |Ψ⟩\,.
```

The exact same form is also assumed in open quantum systems,

```math
i \frac{∂}{∂t} ρ̂ = L(t) ρ̂\,,
```

where ``L`` is the Liouvillian up to a factor of ``i``, see [`liouvillian`](@ref) with `convention=:TDSE`.

## Control Functions

While the generator may have an explicit time dependency, in the context of optimal control the time dependency will usually be implicit in a dependency of the generator on one or more control functions. That is, ``Ĥ(t) = Ĥ(\{ϵ_l(t)\}, t)``. Any object passed to [`propagate`](@ref) as a `generator` must implement the following methods:

* [`QuantumPropagators.Generators.getcontrols`](@ref) — extract the controls ``\{ϵ_l(t)\}`` from ``Ĥ(\{ϵ_l(t)\}, t)``
* [`QuantumPropagators.Generators.evalcontrols`](@ref), [`QuantumPropagators.Generators.evalcontrols!`](@ref) — plug in values for the controls as well as any explicit time dependency to produce a static operator from the time-dependent generator
* [`QuantumPropagators.Generators.substitute_controls`](@ref)  — replace the controls with different (e.g., optimized) controls, producing a new time-dependent generator
* [`QuantumPropagators.Generators.getcontrolderiv`](@ref) — return ``\frac{∂ Ĥ(\{ϵ_{l'}(t)\}, t)}{∂ ϵ_l(t)}``. The result may be `nothing` if ``Ĥ(t)`` does not depend on the particular ``ϵ_l(t)``, a static operator if ``Ĥ(t)`` is linear in the control ``ϵ_l(t)`` and has no explicit time dependency, or a new time-dependent generator to be evaluated via [`evalcontrols`](@ref QuantumPropagators.Generators.evalcontrols) otherwise.

When writing a custom generator type, the [`QuantumPropagators.Generators.check_generator`](@ref) routine should be used to check that all of these methods are implemented.

## Built-in [`hamiltonian`](@ref) and [`liouvillian`](@ref) Generators.

While any generator ``Ĥ(t) = Ĥ(\{ϵ_l(t)\}, t)`` can be implemented via a custom type, the built-in [`hamiltonian`](@ref) and [`liouvillian`](@ref) initialize a [`Generator`](@ref QuantumPropagators.Generators.Generator) object encapsulating a generator of the form

```math
Ĥ(t) = Ĥ_0 + \sum_l a_l(\{ϵ_{l'}(t)\}, t) \, Ĥ_l\,,
```

that is, an (optional) drift term ``Ĥ₀`` and an arbitrary number of control terms that separate into a scalar control amplitude ``a_l(\{ϵ_{l'}(t)\}, t)`` and a static control operator  ``Ĥ_l``.

## Control Amplitudes

In the form represented by a [`Generator`](@ref QuantumPropagators.Generators.Generator), the time-dependence of the control terms is via control amplitudes ``a_l(\{ϵ_{l'}(t)\}, t)``, which may depend on one or more control function ``\{ϵ_{l'}(t)\}``. In most cases, ``a_l(t) ≡ ϵ_l(t)``. Any other dependency of the control amplitudes on the control functions must be implemented via a custom type, for which the following methods must be defined:

* [`QuantumPropagators.Generators.getcontrols`](@ref)
* [`QuantumPropagators.Generators.evalcontrols`](@ref)
* [`QuantumPropagators.Generators.substitute_controls`](@ref)
* [`QuantumPropagators.Generators.getcontrolderiv`](@ref)

These are similar to the equivalent methods for a custom generator. However, [`QuantumPropagators.Generators.evalcontrols`](@ref) for an amplitude returns a number, not an operator. Also, [`QuantumPropagators.Generators.getcontrolderiv`](@ref) returns `0.0` if the control amplitude does not depend on a particular control function, instead of `nothing`.

Any custom amplitude implementation should be checked with [`QuantumPropagators.Generators.check_amplitude`](@ref)


## Operators

The [`QuantumPropagators.Generators.evalcontrols`](@ref) method applied to a generator returns a static operator. For any operator object the 5-argument [`LinearAlgebra.mul!`](https://docs.julialang.org/en/v1/stdlib/LinearAlgebra/#LinearAlgebra.mul!) must be implemented to apply the operator to a state.

For a [`Generator`](@ref) instance obtained from [`hamiltonian`](@ref) or [`liouvillian`](@ref), the resulting operator will be an [`Operator`](@ref QuantumPropagators.Generators.Operator) instance. Any custom type intended as an operator should use [`QuantumPropagators.Generators.check_operator`](@ref) to verify the implementation.

In practice, an operator might also be used as a generator (for a propagation without any time dependency). For any custom operator type that should support this, the implementation should also be checked with [`QuantumPropagators.Generators.check_generator`](@ref).
