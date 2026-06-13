<!--
SPDX-FileCopyrightText: © 2025 Michael Goerz <mail@michaelgoerz.net>

SPDX-License-Identifier: CC-BY-4.0
-->

# Release Notes

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.0.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

Also see the [GitHub Releases](https://github.com/JuliaQuantumControl/QuantumPropagators.jl/releases).

## [Unreleased]

* Added: `ExponentialUtilitiesPropagator`, a propagator based on `expv` from [ExponentialUtilities.jl](https://github.com/SciML/ExponentialUtilities.jl), with Krylov-subspace caching for in-place propagation [[#97]]
* Added: `GuidedAmplitude`, a control amplitude defined relative to an existing reference amplitude [[#111]]
* Added: `size` and `eltype` to the operator interface [[#86], [#98]]
* Added: Trait-based interface introspection with the new `supports_matrix_interface` and `supports_vector_interface` traits for operators [[#99], [#100], [#101]]
* Changed: `supports_inplace` is now a type-based trait instead of acting on values [[#99]]
* Changed: `LockedAmplitude` and `ShapedAmplitude` were refactored so that the shape and the control may independently be vectors or functions, allowing optimized pulses to be substituted into existing amplitudes [[#110]]
* Changed: Discretized controls are now guaranteed to be `Vector{Float64}` [[#110]]
* Changed: The operator interface no longer requires `setindex!`, which is not implementable for a lazy sum [[#100]]
* Changed: `check_state` now requires that `convert(typeof(state), similar(state))` is defined [[#108]]
* Changed: Automatic spectral-range estimation now silently falls back to Arnoldi for operators that are not `Array`, instead of emitting a warning
* Fixed: The fallback for `supports_inplace` on `AbstractArray`/`AbstractVector` was obsolete and inconsistent [[#103]]

## [v0.8.5] — 2025-10-30

* Fixed: The Newton propagator no longer relies on broken Julia internals [[#91], [#93]]

## [v0.8.4] — 2025-10-07

* Added: `CRABFunction` and `VariedFrequencyCRABFunction` parameterized controls, implementing the Chopped Random Basis (CRAB) ansatz, in a new `ParameterizedFunctions` submodule [[#90]]
* Fixed: Removed a type instability in the Newton propagator

## [v0.8.3] — 2024-11-13

* Added: Timing benchmarks for `ExpProp` [[#84]]
* Fixed: Tests now pass on Julia 1.11

## [v0.8.2] — 2024-09-24

* Fixed: The callback is now called before observables, matching the documented behavior [[#83]]

## [v0.8.1] — 2024-09-03

* Removed: The dependency on `QuantumControlBase`; `QuantumPropagators` is now self-contained

## [v0.8.0] — 2024-07-27

* Added: `supports_inplace` to check whether a state or operator type is suitable for in-place operations
* Changed: The in-place versus not-in-place interface is now defined explicitly via the behavior guaranteed by `check_state`, `check_generator`, and `check_operator`, rather than relying on `ismutable`
* Added: Documentation for in-place propagation

## [v0.7.6] — 2024-05-19

* Fixed: Documentation links

## [v0.7.5] — 2024-04-21

* Added: The `t_mid` function is now exposed
* Changed: Improved support for immutable states and operators

## [v0.7.4] — 2024-04-17

* Added: `propagate_sequence` for sequential propagations with different parameters
* Changed: Generalized and extended `get_parameters`
* Changed: The `try`/`catch` in `_make_generator` is now hidden behind the `check` flag

## [v0.7.3] — 2024-01-22

* Changed: Revised `show` for generators and operators

## [v0.7.2] — 2024-01-08

* Changed: Documentation improvements

## [v0.7.1] — 2023-12-19

* Added: `QuantumPropagators.VERSION`
* Fixed: `Op * SVector`

## [v0.7.0] — 2023-12-13

* Added: Support for ODE propagation via OrdinaryDiffEq.jl, as an optional dependency
* Added: Parametrized controls, with a `get_parameters` routine and a `ParameterizedFunction` abstract type
* Added: The spectral range for the Cheby propagator can now be specified manually, via `E_min`/`E_max` and `specrange_method`
* Added: `preserve_start` and `preserve_end` options to `get_tlist_midpoints`
* Added: `check_tlist` interface routine
* Changed: The propagation `method` must now be specified explicitly as a `Module`; the brittle `:auto` heuristic was removed
* Changed: `propagate` was split into lower-level components; `propagate(propagator)` can be called as a low-overhead option to propagate over an entire time grid
* Changed: All propagators now store `timer_data`
* Changed: Interface-check routines now show a truncated backtrace to make problems easier to pinpoint

## [v0.6.1] — 2023-10-04

* Added: Internal profiling via [TimerOutputs.jl](https://github.com/KristofferC/TimerOutputs.jl), with a `timings_enabled` function
* Added: `check_propagator` routine
* Changed: The minimum supported Julia version is now 1.9
* Changed: Refactored `discretize_on_midpoints`

## [v0.6.0] — 2023-05-15

* Added: `Interfaces` submodule with interface-checking routines
* Added: `check` argument to `propagate`
* Changed: Relaxed the requirement of `similar` for states

## [v0.5.0] — 2023-04-04

* Added: Support for time-dependent observables
* Changed: Simplified the public `write_to_storage!` API to `write_to_storage!(storage, i, data)`

## [v0.4.2] — 2023-03-18

* Fixed: The discretization routines now ensure a copy of their input

## [v0.4.1] — 2023-03-15

* Fixed: `NaN` in `flattop` with `t_rise=0` [[#42]]

## [v0.4.0] — 2023-02-15

* Changed: `set_state!` is now private

## [v0.3.1] — 2023-01-25

* Added: `uniform_dt_tolerance` for the Cheby propagator
* Changed: Improved the error message for `evaluate`

## [v0.3.0] — 2022-12-01

* Added: Lazy `Generator` and `Operator` types that apply a time-dependent generator term-by-term to a state instead of summing into a single matrix; significantly faster, especially for sparse matrices
* Added: `hamiltonian` and `liouvillian` functions
* Added: Submodules `Amplitudes`, `Controls`, and `Shapes` (moved from `QuantumControlBase`)
* Added: The concept of "control amplitudes" within the new `Generator` type
* Changed: Renamed methods to use underscores: `initprop` → `init_prop`, `propstep!` → `prop_step!`, `reinitprop!` → `reinit_prop!`, `getcontrols` → `get_controls`, and `getcontrolderiv`/`getcontrolderivs` → `get_control_deriv`/`get_control_derivs`
* Changed: Renamed `eval_controls` → `evaluate` and generalized it to evaluate arbitrary time-dependent objects to static objects; the `vals_dict` for plugging in control values is now an optional keyword argument
* Changed: Renamed `substitute_controls` → `substitute` and generalized it to a wide range of objects
* Changed: `get_control_deriv` now returns a generator, which can be converted into an `Operator` via `evaluate`

## [v0.2.1] — 2022-09-25

* Changed: Refactored to support explicitly time-dependent generators
* Fixed: Dispatch ambiguity

## [v0.2.0] — 2022-09-08

* Added: Initial implementation of the Propagator interface [[#26]]
* Removed: `get_control_parameters` (temporarily)
* Fixed: `axpy!` incompatibility with Julia 1.9

## [v0.1.6] — 2022-03-30

* Fixed: `specrange` now returns `Float64`

## [v0.1.5] — 2022-03-29

* Added: Automatic choice of Cheby for generators with real eigenvalues
* Changed: Minor performance improvements

## [v0.1.4] — 2022-03-23

* Fixed: Spectral radius via diagonalization (`specrad`)

## [v0.1.3] — 2022-03-22

* Added: Option to print an info message in `initpropwrk`

## [v0.1.2] — 2022-02-27

* Fixed: Use of Cheby for backward propagation

## [v0.1.1] — 2022-02-24

* Fixed: Calculation of the spectral range for Cheby
* Fixed: Calculation of expectation values [[#13]]

## [v0.1.0] — 2022-02-06

Initial public release

[Unreleased]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/compare/v0.8.5..HEAD
[v0.8.5]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/releases/tag/v0.8.5
[v0.8.4]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/releases/tag/v0.8.4
[v0.8.3]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/releases/tag/v0.8.3
[v0.8.2]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/releases/tag/v0.8.2
[v0.8.1]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/releases/tag/v0.8.1
[v0.8.0]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/releases/tag/v0.8.0
[v0.7.6]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/releases/tag/v0.7.6
[v0.7.5]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/releases/tag/v0.7.5
[v0.7.4]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/releases/tag/v0.7.4
[v0.7.3]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/releases/tag/v0.7.3
[v0.7.2]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/releases/tag/v0.7.2
[v0.7.1]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/releases/tag/v0.7.1
[v0.7.0]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/releases/tag/v0.7.0
[v0.6.1]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/releases/tag/v0.6.1
[v0.6.0]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/releases/tag/v0.6.0
[v0.5.0]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/releases/tag/v0.5.0
[v0.4.2]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/releases/tag/v0.4.2
[v0.4.1]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/releases/tag/v0.4.1
[v0.4.0]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/releases/tag/v0.4.0
[v0.3.1]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/releases/tag/v0.3.1
[v0.3.0]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/releases/tag/v0.3.0
[v0.2.1]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/releases/tag/v0.2.1
[v0.2.0]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/releases/tag/v0.2.0
[v0.1.6]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/releases/tag/v0.1.6
[v0.1.5]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/releases/tag/v0.1.5
[v0.1.4]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/releases/tag/v0.1.4
[v0.1.3]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/releases/tag/v0.1.3
[v0.1.2]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/releases/tag/v0.1.2
[v0.1.1]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/releases/tag/v0.1.1
[v0.1.0]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/releases/tag/v0.1.0
[#13]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/pull/13
[#26]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/issues/26
[#42]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/issues/42
[#83]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/pull/83
[#84]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/pull/84
[#86]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/issues/86
[#90]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/pull/90
[#91]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/issues/91
[#93]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/pull/93
[#97]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/pull/97
[#98]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/pull/98
[#99]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/pull/99
[#100]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/pull/100
[#101]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/pull/101
[#103]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/pull/103
[#108]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/pull/108
[#110]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/pull/110
[#111]: https://github.com/JuliaQuantumControl/QuantumPropagators.jl/pull/111
