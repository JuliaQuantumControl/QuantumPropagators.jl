```@meta
# SPDX-FileCopyrightText: © 2026 Michael Goerz <mail@michaelgoerz.net>
# SPDX-License-Identifier: MIT OR CC-BY-4.0
```

# [Time-Dependent Controls](@id TimeDependentControls)

```@example control-diagrams
using Main.ControlDiagrams: diagram, legend  # hide
nothing  # hide
```

QuantumPropagators is a framework for simulating dynamics, but its data structures and concepts are rooted in the context of quantum control. Most importantly, we impose a structure on how [dynamical generators](@ref "Dynamical Generators") express their time dependencies as [_amplitudes_](@extref QuantumControl Control-Amplitude) and [_controls_](@extref QuantumControl Control-Function), where a control is a function ``ϵ_l(t)`` that an optimal control method tunes to steer the dynamics towards a particular goal, with the index ``l`` numbering the different controls that the Hamiltonian might depend on.

There are fundamentally two categories of methods:

* Treat ``ϵ_l(t)`` as an arbitrary function that is tuned as a whole. Practically, ``ϵ_l(t)`` is discretized to its values on a time grid (as a ["pulse"](@extref QuantumControl Pulse)). This is what [Krotov's method](@extref Krotov :doc:`index`) and [GRAPE](@extref GRAPE :doc:`index`) do.
* Tune a few scalar [control parameters](@extref QuantumControl Control-Parameters) ``u_n`` inside an analytic ``ϵ_l(\{u_n\}, t)``. This is what gradient-free optimization with, e.g., a Nelder-Mead method does, and what gradient-based methods like GOAT need.

Both are represented in a unified system. A [generator](@ref "Dynamical Generators") (a Hamiltonian or Liouvillian) contains [control amplitudes](@ref ControlAmplitudes), the amplitudes contain controls, and the controls may hold parameters. The diagrams below illustrate that structure and what the QuantumPropagators package provides. This includes [_evaluating_](@ref QuantumPropagators.Controls.evaluate) controls, amplitudes, and generators at a particular point in time and/or by plugging in a particular value for ``ϵ_l(t)`` (independent of what the actual value of the function is at ``t``, so that we can actually _vary_ it). It also includes _discretizing_ [on a time grid](@ref QuantumPropagators.Controls.discretize) or [on the intervals of a time grid](@ref QuantumPropagators.Controls.discretize_on_midpoints). With controls that are defined in terms of analytic parameters, we have functions to recursively extract and collate those parameters in such a way that a control method has direct access to vary the parameters in-place, even if they are deeply nested.

The page [Controls and Derivatives](@extref QuantumControl ControlsAndDerivatives) in the documentation of QuantumControl shows the same diagrams, extended by the additional tools that gradient-based optimization methods require: the derivatives with respect to [a control function as a whole](@extref `QuantumControl.Controls.get_control_deriv`) and with respect to [a particular control parameter](@extref `QuantumControl.Controls.get_parameter_deriv`).


## [Control Functions](@id ControlFunctions)

```@example control-diagrams
legend(:functions; part = :top)  # hide
```

```@example control-diagrams
diagram(:control_functions; part = :top)  # hide
```

The built-in [`Generator`](@ref QuantumPropagators.Generators.Generator) ``Ĥ = Ĥ_0 + ϵ_1(t) Ĥ_1 + a_2(t) Ĥ_2 + a_3(t) Ĥ_3`` in the diagram, as created by [`hamiltonian`](@ref QuantumPropagators.Generators.hamiltonian), has one term of each kind. The control ``ϵ_1`` is its own amplitude. The amplitude ``a_2`` is a [`ShapedAmplitude`](@ref QuantumPropagators.Amplitudes.ShapedAmplitude), a control ``ϵ_2`` times a fixed shape ``S(t)``. The amplitude ``a_3 = ϵ_3^2(t)`` is a custom type that is non-linear in its control.

The terms in this structure are:

* A *control* ``ϵ_l(t)`` is a function of ``t``, or a pulse. The controls of an object are returned by [`get_controls`](@ref QuantumPropagators.Controls.get_controls).
* A *control amplitude* ``a_l(t)`` is the coefficient of the control operator ``Ĥ_l`` in the generator. It can be the control itself, or depend on several controls and on time.
* A *pulse* is a control discretized to the midpoints of the time grid: a `Vector{Float64}` with one value for each of the ``N_T`` intervals of the time grid `tlist`. A pulse is a valid control itself.

The following methods connect the objects in the diagram:

* [`discretize_on_midpoints(ϵ, tlist)`](@ref QuantumPropagators.Controls.discretize_on_midpoints) returns the pulse for the control `ϵ`.
* [`discretize(ϵ, tlist)`](@ref QuantumPropagators.Controls.discretize) returns the values of `ϵ` on the ``N_T + 1`` points of `tlist`. Krotov's method and GRAPE use it to report the optimized controls.
* [`evaluate(obj, tlist, n; vals_dict)`](@ref QuantumPropagators.Controls.evaluate) evaluates `obj` at the midpoint ``t_n`` of the `n`-th interval of `tlist`. For a generator, the result is the operator ``Ĥ_0 + \sum_l a_l(t_n) Ĥ_l``. For an amplitude or a control, it is a `Float64`. A static operator such as ``Ĥ_1`` evaluates to itself.

An entry `ϵ => value` in `vals_dict` replaces the value of the control `ϵ` in the evaluation. This is how the piecewise-constant propagators that Krotov's method and GRAPE use evaluate the generator: with `vals_dict = IdDict(ϵₗ => pₗ[n])`, where `pₗ` is the pulse that the [propagator](@ref QuantumPropagators.AbstractPropagator) stores for the control `ϵₗ`. The optimization changes the pulses, and never the generator.


## [Control Parameters](@id ControlParameters)

A control with parameters is typically an instance of a [`ParameterizedFunction`](@ref QuantumPropagators.Controls.ParameterizedFunction): a [function-like object](@extref Julia Function-like-objects) of `t` with a field `parameters`. See the how-to on [defining a parameterized control](@ref howto_parameterized) for an example. The parameters of any object are returned by [`get_parameters`](@ref QuantumPropagators.Controls.get_parameters), see the how-to on [working with control parameters](@ref howto_parameters). The vectors that store the parameters are [parameter blocks](@ref ParameterBlocks), which the objects in the following diagrams *hold* or *inherit*.

```@example control-diagrams
legend(:parameters; part = :top)  # hide
```

For the built-in [`Generator`](@ref QuantumPropagators.Generators.Generator), only the controls hold parameters:

```@example control-diagrams
diagram(:parameters_generator; part = :top)  # hide
```

Two drives from one source share the carrier frequency ``ω``, but have separate amplitudes ``A_1`` and ``A_2``. Therefore, ``ϵ_1`` and ``ϵ_2`` hold the same block `u = [ω, A₁, A₂]`, and neither is its owner. Each reads only the entries it needs.

The result of `get_parameters(H)` combines the blocks `u` and `v` and stores no values of its own. Its index 2 is `u[2]`, and its index 4 is `v[1]`. Writing to it writes into `u` and `v`.

A custom generator type is not bound to this structure. It can hold blocks in its own fields, and inherit the blocks of the controls that it contains:

```@example control-diagrams
diagram(:parameters_custom; part = :top)  # hide
```

The `TransmonGenerator` in this diagram describes a transmon with a tunable frequency ``ν`` and anharmonicity ``α`` in the drift. The generator `G` holds the block `w` in one of its fields, and it contains a Gaussian control created outside of it. Its own [`get_parameters`](@ref QuantumPropagators.Controls.get_parameters) method combines `w` with the blocks of its control, using [`combine_parameters`](@ref QuantumPropagators.Controls.combine_parameters):

```julia
get_parameters(G::TransmonGenerator) = combine_parameters(G.w, get_parameters(G.control))
```

The order of the blocks is up to that method. This one puts its own block first. The how-tos implement the [`GaussianControl`](@ref howto_parameterized) and the [`TransmonGenerator`](@ref howto_custom_parameters).

The built-in `Generator` cannot hold blocks itself. To add tunable parameters to its drift (or to a dissipator), use [`StaticParameter`](@ref QuantumPropagators.Controls.StaticParameter), a control that is constant in time and holds a block. For the transmon above,

```julia
H = hamiltonian(
    (b⁺b, StaticParameter(w, 1)),
    (b⁺b⁺bb / 2, StaticParameter(w, 2)),
    (b + b⁺, ϵ₄)
)
```

creates a `Generator` with the same parameters as `G`.


### [Parameter Blocks](@id ParameterBlocks)

The diagrams use the following terms:

* A *parameter block* is a vector of parameters that is stored as one object, typically the `parameters` field of a control. It is the unit of identity (`===`) and of sharing: two controls share parameters by holding the same block.
* An object *holds* a block when the block is one of its fields and its own [`get_parameters`](@ref QuantumPropagators.Controls.get_parameters) method returns it directly. A [`ParameterizedFunction`](@ref QuantumPropagators.Controls.ParameterizedFunction) holds exactly one block. A plain function or a pulse holds none. Several objects can hold the same block, and none of them is special: a block usually exists before any of them, created by the user and passed to their constructors.
* An object *inherits* the blocks held by the controls that it contains. An amplitude and the built-in `Generator` only inherit. A custom generator type can do both. Inheriting has no arrow of its own in the diagrams: it is the path through "contains" to a holder.
* *The parameters of an object* are what `get_parameters(object)` returns: an `ArrayPartition` of every block that the object holds or inherits, each once, in order of first occurrence.
* A *parameter* is the pair `(block, i)` of a block and a linear index into it.

The parameters of an object are a [`RecursiveArrayTools.ArrayPartition`](@extref) also if the object has a single block, or none, so that all objects have the same kind of parameters. Each block occupies a contiguous range of indices in it, like the blocks of a `BlockArrays.BlockVector`, but the blocks remain separate arrays in memory. The `ArrayPartition` stores no values of its own, so writing to it writes into the blocks. It is neither a block nor a holder. As long as an object holds and contains the same blocks, `get_parameters` returns identical (`===`) results.

The field `x` of the parameters is the tuple of blocks, in the order of their index ranges, so that a linear index into the parameters translates into a pair `(block, i)`.

A block must meet four requirements, which are verified by [`check_parameterized`](@ref QuantumPropagators.Interfaces.check_parameterized) (the second one only partially):

1. It is an `AbstractVector{Float64}` with one-based indexing.
2. It can be mutated in place, and every holder reads its values when it is evaluated, so that writing to it changes what [`evaluate`](@ref QuantumPropagators.Controls.evaluate) returns.
3. It is not an `ArrayPartition`. That type is reserved for the parameters of an object.
4. Any two blocks are either `===` or share no memory. Controls that share some of their parameters hold one block that contains all of them, like the block `u` in the [first diagram for control parameters](@ref ControlParameters).

A violation of the fourth requirement does not crash and does not corrupt memory, but gradient-based optimization is silently wrong. If two blocks overlap, the same parameter appears twice in the parameters of the problem. Whatever an optimizer writes to the first of the two entries is overwritten when it writes to the second, and each of the two gradient entries contains only part of the true derivative. A gradient-free optimizer still runs, but searches an extra dimension and reports misleading parameter values.
