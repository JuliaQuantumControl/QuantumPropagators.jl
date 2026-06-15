<!--
SPDX-FileCopyrightText: © 2021 Michael Goerz <mail@michaelgoerz.net>

SPDX-License-Identifier: MIT OR CC-BY-4.0
-->

# QuantumPropagators

[![version](https://juliahub.com/docs/General/QuantumPropagators/stable/version.svg)](https://juliahub.com/ui/Packages/General/QuantumPropagators)
[![Stable](https://img.shields.io/badge/docs-stable-blue.svg)](https://juliaquantumcontrol.github.io/QuantumPropagators.jl/stable)
[![Dev](https://img.shields.io/badge/docs-dev-blue.svg)](https://juliaquantumcontrol.github.io/QuantumPropagators.jl/dev)
[![Build Status](https://github.com/JuliaQuantumControl/QuantumPropagators.jl/workflows/CI/badge.svg)](https://github.com/JuliaQuantumControl/QuantumPropagators.jl/actions)
[![REUSE](https://github.com/JuliaQuantumControl/QuantumPropagators.jl/actions/workflows/REUSE.yml/badge.svg)](https://github.com/JuliaQuantumControl/QuantumPropagators.jl/actions/workflows/REUSE.yml)
[![Coverage](https://codecov.io/gh/JuliaQuantumControl/QuantumPropagators.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/JuliaQuantumControl/QuantumPropagators.jl)

Methods for simulating the time dynamics of a quantum system for packages with the [JuliaQuantumControl][] organization.

## Installation

For standalone use, the `QuantumPropagators` package can be installed with [Pkg][] as

~~~
pkg> add QuantumPropagators
~~~

Otherwise, the package should be installed as part of [`QuantumControl.jl`][QuantumControl]:

~~~
pkg> add QuantumControl
~~~

For development usage, see the [organization development notes](https://github.com/JuliaQuantumControl#development).

## Documentation

The documentation of `QuantumPropagators.jl` is available at <https://juliaquantumcontrol.github.io/QuantumPropagators.jl>.

For a broader perspective on how time propagation fits in the context of quantum control theory, see the [documentation of the `QuantumControl.jl` meta-package](https://juliaquantumcontrol.github.io/QuantumControl.jl/).

## License

This project is [REUSE-compliant](https://reuse.software). The entire package can be distributed under the [MIT License](LICENSE): every file is available under MIT. Files that are not source code are additionally dual-licensed for reuse in other contexts — documentation (including the `README` and other prose) under [Creative Commons Attribution (`CC-BY-4.0`)](https://creativecommons.org/licenses/by/4.0/deed.en), and trivial configuration or tooling files as public domain ([`CC0-1.0`](https://creativecommons.org/publicdomain/zero/1.0/deed.en)). The authoritative per-file copyright and licensing information is provided by SPDX tags in each file (and `REUSE.toml` where a file cannot carry a header). It can be verified using the [`reuse` tool](https://github.com/fsfe/reuse-tool?tab=readme-ov-file#reuse), e.g., by running `reuse lint`.

[QuantumControl]: https://github.com/JuliaQuantumControl/QuantumControl.jl#readme
[JuliaQuantumControl]: https://github.com/JuliaQuantumControl
[Pkg]: https://pkgdocs.julialang.org/v1/
