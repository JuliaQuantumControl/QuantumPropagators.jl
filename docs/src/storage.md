# Expectation Values

The [`propagate`](@ref) routine allows the storage of data for every point of the time grid.  This is done by passing it a `storage` object created with [`QuantumPropagators.Storage.init_storage`](@ref), or simply `storage=true` in order to create the appropriate storage automatically.

By default, the `storage` will be used to store the propagated states at each point in time. However, more generally, what goes into `storage` can be customized via the `observables` parameter of [`propagate`](@ref). This allows to obtain arbitrary expectation values from the propagation. In short, the following options are available:

* Pass `observables` as a tuple of functions that each take the propagated state `Ψ` and return some expectation value.
* For time-dependent expectation values, the functions may take positional parameters `tlist` and `i` in addition to `Ψ`, that is, `(Ψ, tlist, i) -> value`, for an expectation value at time `tlist[i]`.
* Pass `observables` as a tuple of matrices for which to calculate expectation values via the three-argument `dot` function. This assumes the propagated state to be a vector.
* For maximum flexibility, pass `observables` as a tuple of arbitrary objects for which a custom [`map_observable`](@ref QuantumPropagators.Storage.map_observable) method has been defined, or pass `observables` as a single object for which a custom [`map_observables`](@ref QuantumPropagators.Storage.map_observables) method has been defined.


### Details

After each propagation step, with a propagated state at time slot `i`,

1. [`data = QuantumPropagators.Storage.map_observables(observables, tlist, i, state)`](@ref QuantumPropagators.Storage.map_observables) generates `data` from the propagated state
2. [`QuantumPropagators.Storage.write_to_storage!(storage, i, data)`](@ref QuantumPropagators.Storage.write_to_storage!) places that `data` into `storage` for time slot `i`

After [`propagate`](@ref) returns, the [`QuantumPropagators.Storage.get_from_storage!`](@ref) routine can be used to extract data from any time slot. This interface hides the internal memory organization of `storage`, which is set up by [`init_storage`](@ref QuantumPropagators.Storage.init_storage) based on the type of `state` and the given `observables`. This system can be extended with multiple dispatch, allowing to optimize the `storage` for custom data types. Obviously, [`init_storage`](@ref QuantumPropagators.Storage.init_storage), [`map_observables`](@ref QuantumPropagators.Storage.map_observables), [`write_to_storage!`](@ref QuantumPropagators.Storage.write_to_storage!), and [`get_from_storage!`](@ref QuantumPropagators.Storage.get_from_storage!) must all be consistent.

The default implementation of these routine uses either a standard Vector or a Matrix as `storage`.

Roughly speaking, when storing states, if the state of some arbitrary type,
the storage will be a Vector where the i'th entry points to a copy of the propagated state at the i'th time slot. If the state is a Vector, the storage will be a Matrix containing the state for the i'th time slot in the i'th column.

When a tuple of `observables` is passed to [`propagate`](@ref), if [`map_observables`](@ref QuantumPropagators.Storage.map_observables) returns data of the same type for each observable, `storage` will be a Matrix containing the values from the different observables for the i'th time slot in the i'th column. This would be typical for the storage of expectation values, e.g. with

~~~julia
observables=(state->dot(state, Ô₁, state), state->dot(state, Ô₂, state))
~~~

where `Ô₁`, `Ô₂` are two Hermitian operators. In this case, `storage` would be
a `2 × nt` `Float64` array. Calling `get_from_storage!(data, storage, i)` would be
equivalent to `copyto!(data, storage[:,i])` and extract the i'th column of `storage`, i.e. the
vector `[⟨Ô₁⟩, ⟨Ô₂A]⟩]` at time slot `i`. Alternatively, `storage[1,:]` would return the values of `⟨Ô₁⟩` over time. This would be useful for plotting, and illustrates the benefits of using a Matrix as `storage`.

Usually, the `observables` should be functions acting on the `state`, but [`map_observables`](@ref QuantumPropagators.Storage.map_observables) can be extended to other types of observables as well. For example, the situation were `state` is a vector and the `observables` are matrices is also supported; if  `Ô₁`, `Ô₂` are matrices,
~~~julia
observables=(Ô₁, Ô₂)
~~~
would have the same result of storing the expectation values of those two operators.

The `observables` are not required to yield real values: the term "observable"
is used very loosely here. We could directly calculate, e.g., the complex
amplitude α of a coherent state in quantum optics, or the number of levels with
non-zero population (as an integer).

It is possible to have time-dependent "observables", for example, to store "lab-frame" states from a propagation in a (time-dependent) [rotating frame](https://en.wikipedia.org/wiki/Rotating-wave_approximation). This is supported by default by passing a function with the three arguments `state`, `tlist`, `i` as an observable, where `state` is defined to be at time `tlist[i]`.

If there are multiple observables that return data of different types, by default `storage` will be a Vector that contains tuples with the result for each observable. For example, with

~~~julia
observables=(state->dot(state, Ô₁, state), state->count_poplevels(state))
~~~

where `count_poplevels` is a function that counts the number of levels with
non-zero population, the resulting `storage` would be a `Vector{Tuple{Float64, Int64}}`.

If there is a single observable that yields a vector, that vector is stored in the i'th column of a `storage` matrix. This is in fact what happens when storing the propagated states (`observables=(Ψ->copy(Ψ), )`) if `Ψ` is a Vector, but there are other use cases, such as calculating the population in all levels in one go, with
`observables=(Ψ -> abs.(Ψ).^2, )`.

If there is a single variable that yields a non-vector object, `storage` will be a Vector where the i'th entry points to the object. This is in fact what happens by default when storing states  that are e.g. instances of
`QuantumOptics.Ket`. In such a case, it might be advisable to add new methods for [`QuantumPropagators.Storage.init_storage`](@ref) and [`QuantumPropagators.Storage.write_to_storage!`](@ref) that implement a more efficient in-place storage.
