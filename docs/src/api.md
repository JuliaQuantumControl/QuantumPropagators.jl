# API

## [`QuantumPropagators`](@id QuantumPropagatorsAPI)

The `QuantumPropagators.jl` package exports the following high-level routines:

* [`initpropwrk`](@ref) — Initialize a work space for propagation
* [`propagate`](@ref) — Propagate a state over an entire time grid.
* [`propstep`](@ref) — Perform a single propagation step
* [`propstep!`](@ref) — Perform a single propagation step in-place.

These delegate to routines implementing the various propagation method defined in the respective [sub-modules](@ref QuantumPropagatorsChebyAPI).

Furthermore, the following routines from [the `Storage` sub-module](@ref QuantumPropagtorsStorageAPI) are re-exported:

* [`init_storage`](@ref) – Create a storage array for propagated states or expectation values
* [`write_to_storage!`](@ref) – Place data into storage
* [`get_from_storage!`](@ref) – Obtain data from storage (in-place)


### Reference

```@autodocs
Modules = [QuantumPropagators]
```

## [`QuantumPropagators.Storage`](@id QuantumPropagtorsStorageAPI)

The following routine allow manage and extend storage arrays `storage` parameter in [`propagate`](@ref). See [Storage of states or expectation values](@ref) for more details.
Note that the high-level routines [`init_storage`](@ref), [`write_to_storage!`](@ref), and [`get_from_storage!`](@ref) are also re-exported in the top-level [QuantumPropagators](@ref QuantumPropagatorsAPI) directly.

```@autodocs
Modules = [QuantumPropagators.Storage]
```

## [`QuantumPropagators.Cheby`](@id QuantumPropagatorsChebyAPI)

The following routines implement time propagation via expansion in Chebychev
polynomials.

```@autodocs
Modules = [QuantumPropagators.Cheby]
```

## [`QuantumPropagators.Newton`](@id QuantumPropagtorsNewtonAPI)

The following routines implement time propagation via expansion in Newton
polynomials using a restarted Arnoldi scheme to determine evaluation points.

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

```@autodocs
Modules = [QuantumPropagators.ExpProp]
Public = true
Private = false
```


## [`QuantumPropagators.SpectralRange`](@id QuantumPropagtorsSpectralRangeAPI)

[Chebychev propagation](@ref QuantumPropagatorsChebyAPI) relies on estimating the spectral range of the Hamiltonian, which in turn may be done via [Arnoldi iteration](@ref QuantumPropagtorsArnoldiAPI).

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
