```@meta
CurrentModule = QuantumPropagators
```

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



## Contents

```@contents
Pages=[
 "overview.md",
 "generators.md",
 "methods.md",
 "storage.md",
 "examples/index.md",
 "howto.md",
 "benchmarks.md",
 "api/quantumpropagators.md",
]
Depth = 2
```

## History

See the [Releases](https://github.com/JuliaQuantumControl/QuantumPropagators.jl/releases) on Github.
