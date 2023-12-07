# Profiling Howto

To choose an appropriate propagation method and parameters for a given problem, it is essential to benchmark and profile the propagation.

Consider the following simple example:

```@example profiling
using QuantumControlTestUtils.RandomObjects: random_dynamic_generator, random_state_vector

tlist = collect(range(0, step=1.0, length=101));
N = 200;  # size of Hilbert space
H = random_dynamic_generator(N, tlist);
Ψ₀ = random_state_vector(N);
nothing  # hide
```


## BenchmarkTools

The first line of defense is the use of [`BenchmarkTools`](https://juliaci.github.io/BenchmarkTools.jl/stable/). The `@benchmark` macro allows to generate statistics on how long a call to [`propagate`](@ref) takes.

### Chebychev propagation

For example, we can time the propagation with the Chebychev method:


```@example profiling
using BenchmarkTools
using QuantumPropagators
using QuantumPropagators: Cheby

@benchmark propagate($Ψ₀, $H, $tlist; method=Cheby) samples=10
```

### Newton propagation

Or, the same propagation with the Newton method:

```@example profiling
using QuantumPropagators: Newton

@benchmark propagate($Ψ₀, $H, $tlist; method=Newton) samples=10
```

The result in this case illustrates the significant advantage of the Chebychev method for systems of moderate to small size and unitary dynamics.

When using custom data structures for the dynamical generators or states, `@benchmark` should also be used to optimize lower-level operations as much as possible, e.g. the application of the Hamiltonian to the state.


## TimerOutputs

A lot more insight into the internals of a [`propagate`](@ref) call can be obtained by collecting timing data. This functionality is integrated in `QuantumPropagators` and uses the [`TimerOutputs`](https://github.com/KristofferC/TimerOutputs.jl#readme) package internally.

### Enabling the collection of timing data

To enable collecting internal timing data, call [`QuantumPropagators.enable_timings`](@ref):

```@example profiling
QuantumPropagators.enable_timings()
nothing # hide
```

The status of the data collection can be verified with [`QuantumPropagators.timings_enabled`](@ref).


### Chebychev propagation

Since the call to [`QuantumPropagators.enable_timings`](@ref) invalidates existing compiled code, and to avoid the compilation overhead showing up in the timing data, we call [`propagate`](@ref) once to ensure compilation:

```@example profiling
propagate(Ψ₀, H, tlist; method=Cheby);
nothing  # hide
```

In any subsequent propagation, we could access the timing data in a `callback` to [`propagate`](@ref):

```@example profiling
function show_timing_data(propagator, args...)
    if propagator.t == tlist[end]
        show(propagator.timing_data, compact=true)
    end
end

propagate(Ψ₀, H, tlist; method=:cheby, callback=show_timing_data);
nothing  # hide
```

See the [`TimerOutputs` documentation](https://github.com/KristofferC/TimerOutputs.jl#settings-for-printing) for details on how to print the `timing_data`.

Alternatively, without a `callback`:

```@example profiling
propagator = init_prop(Ψ₀, H, tlist; method=Cheby)
for step ∈ 1:(length(tlist)-1)
    prop_step!(propagator)
end
show(propagator.timing_data, compact=true)
```

The reported runtimes here are less important than the number of function calls and the runtime percentages. In this case, the timing data shows that the propagation is dominated by the matrix-vector products (applying the Hamiltonian to the state), as it should. The percentage would go to 100% for larger Hilbert spaces.

### Newton propagation

For the Newton method:

```@example profiling
propagate(Ψ₀, H, tlist; method=Newton);   # recompilation
propagate(Ψ₀, H, tlist; method=Newton, callback=show_timing_data);
nothing  # hide
```

We see here that the Newton propagation requires more matrix-vector products (2000 compared to 1200 for Chebychev), partly because the Newton propagator is "chunked" to `m_max` applications in each "restart" (10 by default, with 2 restarts required to reach machine precision in this case). Moreover, there is significant overhead beyond just matrix-vector multiplication, which will disappear only for significantly larger Hilbert spaces.

### Disabling the collection of timing data

There there is a [small overhead](https://github.com/KristofferC/TimerOutputs.jl#overhead) associated with collecting the timing data, it should not be enabled "in production". To [`QuantumPropagators.disable_timings`](@ref) function undoes the previous [`QuantumPropagators.enable_timings`](@ref):

```@example profiling
QuantumPropagators.disable_timings()
nothing # hide
```

This again will trigger recompilation of any method that was collecting timing data, removing the associated overhead.
