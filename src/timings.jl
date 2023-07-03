using TimerOutputs: enable_debug_timings, disable_debug_timings

"""Enable the collection of `TimerOutputs` data.

```julia
QuantumPropagators.enable_timings()
```

enables certain portions of the package to collect
[`TimerOutputs`](https://github.com/KristofferC/TimerOutputs.jl) internally.
This aids in profiling and benchmarking propagation methods.

Specifically, after `enable_timings()`, for any [`ChebyPropagator`](@ref)
or [`NewtonPropagator`](@ref), timing data will become available in
`propagator.wrk.timing_data` (as a
[`TimerOutput`](https://github.com/KristofferC/TimerOutputs.jl#usage)
instance). This data is reset when the propagator is re-instantiated with
[`init_prop`](@ref) or re-initialized with [`reinit_prop!`](@ref). This makes
the data local to any call of [`propagate`](@ref).

Note that `enable_timings()` triggers recompilation, so
[`propagate`](@ref) should be called at least twice to avoid
compilation overhead in the timing data. There is still a [small
overhead](https://github.com/KristofferC/TimerOutputs.jl#overhead) for
collecting the timing data.

The collection of timing data can be disabled again
with [`disable_timings`](@ref).

Returns [`QuantumPropagators.timings_enabled()`](@ref timings_enabled), i.e.,
`true` if successful.
"""
function enable_timings()
    enable_debug_timings(@__MODULE__)
    enable_debug_timings(Cheby)
    enable_debug_timings(Newton)
    enable_debug_timings(Arnoldi)
    return timings_enabled()
end


"""Check whether the collection of `TimerOutputs` data is active.

```julia
QuantumPropagators.timings_enabled()
```

returns `true` if [`QuantumPropagators.enable_timings()`](@ref
enable_timings) was called, and `false` otherwise or after
[`QuantumPropagators.disable_timings()`](@ref disable_timings).
"""
function timings_enabled()
    enabled = @eval getfield(@__MODULE__, :timeit_debug_enabled)()
    enabled &= @eval getfield(Cheby, :timeit_debug_enabled)()
    enabled &= @eval getfield(Newton, :timeit_debug_enabled)()
    enabled &= @eval getfield(Arnoldi, :timeit_debug_enabled)()
    return enabled
end


"""Disable the collection of `TimerOutputs` data.

```julia
QuantumPropagators.disable_timings()
```

disables the collection of timing data previously enabled with
[`enable_timings`](@ref). This triggers recompilation to completely remove
profiling from the code. That is, there is [zero
cost](https://github.com/KristofferC/TimerOutputs.jl#overhead) when
the collection of timing data is disabled.

Returns [`QuantumPropagators.timings_enabled()`](@ref timings_enabled), i.e.,
`false` if successful.
"""
function disable_timings()
    disable_debug_timings(@__MODULE__)
    disable_debug_timings(Cheby)
    disable_debug_timings(Newton)
    disable_debug_timings(Arnoldi)
    return timings_enabled()
end
