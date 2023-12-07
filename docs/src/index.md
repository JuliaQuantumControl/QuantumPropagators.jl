# QuantumPropagators

```@eval
using Markdown
using Pkg

VERSION = Pkg.dependencies()[Base.UUID("7bf12567-5742-4b91-a078-644e72a65fc1")].version

github_badge = "[![Github](https://img.shields.io/badge/JuliaQuantumControl-QuantumPropagators.jl-blue.svg?logo=github)](https://github.com/JuliaQuantumControl/QuantumPropagators.jl)"

version_badge = "![v$VERSION](https://img.shields.io/badge/version-v$VERSION-green.svg)"

Markdown.parse("$github_badge $version_badge")
```

The [QuantumPropagators](https://github.com/JuliaQuantumControl/QuantumPropagators.jl) package implements methods for simulating the time dynamics of a quantum system. It is the numerical backend for all packages implementing methods of optimal control within the [JuliaQuantumControl](https://github.com/JuliaQuantumControl) organization.

The purpose is this package is to:

* Implement the piecewise-constant propagation methods that have been traditionally used in quantum control, based on an expansion of the time evolution operator into Chebychev polynomials (closed systems) or Newton polynomials (open systems),
* Delegate to third-party packages that implement state-of-the-art methods for time-continuous dynamics ([DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/)).
* Define data structures for "dynamical generators" (time-dependent Hamiltonians or Liouvillians) that are relevant to quantum control. Specifically, this includes a well-defined notion of "control functions" and "control parameters" inside the Hamiltonian.
* Define a  general interface for a time propagators that is suitable for all aspects of quantum control. This includes the ability to efficiently tune control parameters or control functions while the propagation is running (for "sequential updates"), and the ability to propagate both forwards and backwards in time.

The focus on quantum control separates this package from other solutions for just simulating the time dynamics of a quantum system, e.g. within [QuantumOptics.jl](https://qojulia.org) or directly with [DifferentialEquations.jl](https://docs.sciml.ai/DiffEqDocs/stable/). Also, while `QuantumPropagators` defines a general interface for dynamic generators, it is otherwise agnostic about the data structures to represent Hamiltonians or quantum states. This allows for the use of modeling frameworks like [QuantumOptics.jl](https://qojulia.org) or of custom problem-specific data structures.


## Contents

```@contents
Pages=[
 "overview.md",
 "generators.md",
 "methods.md",
 "storage.md",
 "howto.md",
 "benchmarks.md",
 "api/quantumpropagators.md",
]
Depth = 2
```


## History

See the [Releases](https://github.com/JuliaQuantumControl/QuantumPropagators.jl/releases) on Github.
