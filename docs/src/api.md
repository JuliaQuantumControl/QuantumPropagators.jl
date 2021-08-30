# API

## High level routines

The `QuantumPropagators.jl` package provides the following high-level routines:

* [`initpropwrk`](@ref) — Initialize a work space for propagation
* [`propagate`](@ref) — Propagate a state over an entire time grid.
* [`propstep!`](@ref) — Perform a single propagation step in-place.

These delegate the routines implementing the various propagation methods,
detailed in the reference sections below.

See the [Index](@ref) for the full list of routines.

### Reference

```@autodocs
Modules = [QuantumPropagators]
Pages = ["propagate.jl"]
```

## Storage Reference

The following routine allow manage and extend storage arrays `storage` parameter in [`propagate`](@ref). See [Storage of states or expectation values](@ref) for more details.

```@autodocs
Modules = [QuantumPropagators]
Pages = ["storage.jl"]
```

## Chebychev reference

The following routines implement time propagation via expansion in Chebychev
polynomials.

```@autodocs
Modules = [QuantumPropagators]
Pages = ["cheby.jl"]
```

## Newton reference

The following routines implement time propagation via expansion in Newton
polynomials using a restarted Arnoldi scheme to determine evaluation points.

### Public members

```@autodocs
Modules = [QuantumPropagators]
Pages = ["newton.jl"]
Public = true
Private = false
```

### Private members

```@autodocs
Modules = [QuantumPropagators]
Pages = ["newton.jl"]
Public = false
Private = true
```

## Matrix exponentiation reference

The following routines implement time propagation via explicit exponentiation
of the dynamical generator.

```@autodocs
Modules = [QuantumPropagators]
Pages = ["expprop.jl"]
Public = true
Private = false
```


## Spectral Radius and Arnoldi

Chebychev propagation relies on estimating the spectral radius of the
Hamiltonian, which in turn may be done via Arnoldi iteration.

### Public members

```@autodocs
Modules = [QuantumPropagators]
Pages = ["specrad.jl", "arnoldi.jl"]
Public = true
Private = false
```

### Private members

```@autodocs
Modules = [QuantumPropagators]
Pages = ["specrad.jl", "arnoldi.jl"]
Public = false
Private = true
```


## Index

```@index
```
