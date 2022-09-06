# API

## [`QuantumPropagators`](@id QuantumPropagatorsAPI)

The highest-level API `QuantumPropagators.jl` package consists solely of a single function:

* [`propagate`](@ref) — Propagate a quantum state over an entire time grid

At a slightly lower level, propagation of quantum states in encapsulated by [The Propagator Interface](@ref):

* [`initprop`](@ref) — Initialize a `propagator` object, which is of some concrete (method-dependent) sub-type of [`AbstractPropagator`](@ref QuantumPropagators.AbstractPropagator)
* [`reinitprop!`](@ref) — Re-initialize the `propagator`
* [`propstep!`](@ref) — Advance the `propagator`  by a single time step forward or backward
* [`set_state!`](@ref) — Mutate the current quantum `state` of the `propagator`.
* [`set_t!`](@ref QuantumPropagators.set_t!) — Mutate the current time of the `propagator` (**not exported**)


To support the `storage` feature of [`propagate`](@ref), the following routines from [the `Storage` sub-module](@ref QuantumPropagtorsStorageAPI) are re-exported in `QuantumPropagators`:

* [`init_storage`](@ref) – Create a storage array for propagated states or expectation values
* [`write_to_storage!`](@ref) – Place data into storage
* [`get_from_storage!`](@ref) – Obtain data from storage (in-place)


The above list of methods constitutes the entire *public* interface of `QuantumPropagators`. At the lowest level, further functionality is provided by sub-modules like [`QuantumPropagators.Cheby`](@ref QuantumPropagatorsChebyAPI), which defines a standalone API specifically for the Chebychev propagation method.

The full list of sub-modules is

* [`QuantumPropagators.Controls`](@ref QuantumPropagtorsControsAPI)
* [`QuantumPropagators.Storage`](@ref QuantumPropagtorsStorageAPI)
* [`QuantumPropagators.Cheby`](@ref QuantumPropagatorsChebyAPI)
* [`QuantumPropagators.Newton`](@ref QuantumPropagtorsNewtonAPI)
* [`QuantumPropagators.ExpProp`](@ref QuantumPropagtorsExpPropAPI)
* [`QuantumPropagators.SpectralRange`](@ref QuantumPropagtorsSpectralRangeAPI)
* [`QuantumPropagators.Arnoldi`](@ref QuantumPropagtorsArnoldiAPI)


### Reference

#### Public members
```@autodocs
Modules = [QuantumPropagators]
Public = true
Private = false
```

#### Private members

```@autodocs
Modules = [QuantumPropagators]
Public = false
Private = true
```

## [`QuantumPropagators.Controls`](@id QuantumPropagtorsControsAPI)

The following routines define how controls are extracted from the time-dependent Hamiltonians or Liouvillians

```@index
Modules = [QuantumPropagators.Controls]
```

### Public members

```@autodocs
Modules = [QuantumPropagators.Controls]
Public = true
Private = false
```

## [`QuantumPropagators.Storage`](@id QuantumPropagtorsStorageAPI)

The following routines allow to manage and extend storage arrays `storage` parameter in [`propagate`](@ref). See [Storage of states or expectation values](@ref) for more details.
Note that the high-level routines [`init_storage`](@ref), [`write_to_storage!`](@ref), and [`get_from_storage!`](@ref) are also re-exported in the top-level [QuantumPropagators](@ref QuantumPropagatorsAPI) directly.

```@index
Modules = [QuantumPropagators.Storage]
```

### Public members

```@autodocs
Modules = [QuantumPropagators.Storage]
Public = true
Private = false
```

## [`QuantumPropagators.Cheby`](@id QuantumPropagatorsChebyAPI)

The following routines implement time propagation via expansion in Chebychev
polynomials.

```@index
Modules = [QuantumPropagators.Cheby]
```

### Public members

```@autodocs
Modules = [QuantumPropagators.Cheby]
Public = true
Private = false
```

## [`QuantumPropagators.Newton`](@id QuantumPropagtorsNewtonAPI)

The following routines implement time propagation via expansion in Newton
polynomials using a restarted Arnoldi scheme to determine evaluation points.

```@index
Modules = [QuantumPropagators.Newton]
```

### Public members

```@autodocs
Modules = [QuantumPropagators.Newton]
Public = true
Private = false
```

### Private members

```@autodocs
Modules = [QuantumPropagators.Newton]
Public = false
Private = true
```

## [`QuantumPropagators.ExpProp`](@id QuantumPropagtorsExpPropAPI)

The following routines implement time propagation via explicit exponentiation of the dynamical generator.

```@index
Modules = [QuantumPropagators.ExpProp]
```

### Public members

```@autodocs
Modules = [QuantumPropagators.ExpProp]
Public = true
Private = false
```


## [`QuantumPropagators.SpectralRange`](@id QuantumPropagtorsSpectralRangeAPI)

[Chebychev propagation](@ref QuantumPropagatorsChebyAPI) relies on estimating the spectral range of the Hamiltonian, which in turn may be done via [Arnoldi iteration](@ref QuantumPropagtorsArnoldiAPI).

```@index
Modules = [QuantumPropagators.SpectralRange]
```

### Public members

```@autodocs
Modules = [QuantumPropagators.SpectralRange]
Public = true
Private = false
```

### Private members

```@autodocs
Modules = [QuantumPropagators.SpectralRange]
Public = false
Private = true
```

## [`QuantumPropagators.Arnoldi`](@id QuantumPropagtorsArnoldiAPI)

[Arnoldi iteration](https://en.wikipedia.org/wiki/Arnoldi_iteration) is an approximate method to find extremal eigenvalues of a dynamical generator by projecting it into a [Krylov subspace](https://en.wikipedia.org/wiki/Krylov_subspace). It is used to estimate spectral ranges for [Chebychev Propagation](@ref QuantumPropagatorsChebyAPI) and to find appropriate Leja points to support [Newton Propagation](@ref QuantumPropagtorsNewtonAPI)

```@index
Modules = [QuantumPropagators.Arnoldi]
```


### Private members

```@raw html
<!--There are no public routines in Arnoldi-->
```

```@autodocs
Modules = [QuantumPropagators.Arnoldi]
Public = false
Private = true
```


## Index

```@index
```
