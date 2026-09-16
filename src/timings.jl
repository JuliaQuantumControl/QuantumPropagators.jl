# SPDX-FileCopyrightText: © 2021 Michael Goerz <mail@michaelgoerz.net>
#
# SPDX-License-Identifier: MIT

using TimerOutputs: enable_debug_timings, disable_debug_timings

# All modules with `@timeit_debug` sections, including loaded extensions
function _timing_modules()
    modules = Module[@__MODULE__, Cheby, ExpProp, Newton, Arnoldi]
    for name in (:QuantumPropagatorsODEExt, :QuantumPropagatorsExponentialUtilitiesExt)
        extmod = Base.get_extension(@__MODULE__, name)
        isnothing(extmod) || push!(modules, extmod)
    end
    return modules
end

"""Enable the collection of `TimerOutputs` data.

```julia
QuantumPropagators.enable_timings()
```

enables certain portions of the package to collect
[`TimerOutputs`](@extref TimerOutputs readme) internally.
This aids in profiling and benchmarking propagation methods.

Specifically, after `enable_timings()`, for any [`ChebyPropagator`](@ref)
or [`NewtonPropagator`](@ref), timing data will become available in
`propagator.wrk.timing_data` (as a [`TimerOutput`](@extref TimerOutputs usage)
instance). This data is reset when the propagator is re-instantiated with
[`init_prop`](@ref) or re-initialized with [`reinit_prop!`](@ref). This makes
the data local to any call of [`propagate`](@ref).

Note that `enable_timings()` triggers recompilation, so
[`propagate`](@ref) should be called at least twice to avoid
compilation overhead in the timing data. There is still a [small
overhead](@extref TimerOutputs overhead) for collecting the timing data.

This includes the propagators defined in package extensions, e.g., the
`ExponentialUtilitiesPropagator`, if the extension is loaded at the time
`enable_timings()` is called. An extension that is loaded later requires
calling `enable_timings()` again.

The collection of timing data can be disabled again
with [`disable_timings`](@ref).

Returns [`QuantumPropagators.timings_enabled()`](@ref timings_enabled), i.e.,
`true` if successful.
"""
function enable_timings()
    foreach(enable_debug_timings, _timing_modules())
    return timings_enabled()
end


"""Check whether the collection of `TimerOutputs` data is active.

```julia
QuantumPropagators.timings_enabled()
```

returns `true` if [`QuantumPropagators.enable_timings()`](@ref
enable_timings) was called, and `false` otherwise or after
[`QuantumPropagators.disable_timings()`](@ref disable_timings). It also
returns `false` if any loaded package extension does not collect timing data,
e.g., because the extension was loaded after `enable_timings()`.
"""
function timings_enabled()
    # `enable_debug_timings` redefines `timeit_debug_enabled` in each module,
    # so the call must run in the latest world age
    return all(
        mod -> Base.invokelatest(getfield(mod, :timeit_debug_enabled)),
        _timing_modules()
    )
end


"""Disable the collection of `TimerOutputs` data.

```julia
QuantumPropagators.disable_timings()
```

disables the collection of timing data previously enabled with
[`enable_timings`](@ref). This triggers recompilation to completely remove
profiling from the code. That is, there is
[zero cost](@extref TimerOutputs overhead) when
the collection of timing data is disabled.

Returns [`QuantumPropagators.timings_enabled()`](@ref timings_enabled), i.e.,
`false` if successful.
"""
function disable_timings()
    foreach(disable_debug_timings, _timing_modules())
    return timings_enabled()
end
