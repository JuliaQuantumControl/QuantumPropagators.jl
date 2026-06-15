<!--
SPDX-FileCopyrightText: © 2021 Michael Goerz <mail@michaelgoerz.net>

SPDX-License-Identifier: MIT OR CC0-1.0
-->

# Contributing

Everyone is welcome to contribute! You can contribute simply by opening issues to report bugs or request features.

Development of all packages in the [JuliaQuantumControl] organization follows shared contributing guidelines. Please refer to these for instructions on reporting issues, making pull requests, running the tests, building the documentation, code formatting, the changelog, and the release process:

<https://github.com/JuliaQuantumControl/.github/blob/master/CONTRIBUTING.md>


## Licensing

This package is [REUSE-compliant](https://reuse.software): every file declares its copyright and license via SPDX tags, which CI enforces (verify locally with `make reuse`). When adding a new file, include an SPDX header (or a sidecar `*.license` file when the format has no comment syntax) with the copyright `© <year> Michael Goerz <mail@michaelgoerz.net>`, and choose the license by category: source code (`src/`, `ext/`, `test/`) → `MIT`; prose/documentation (`README`, `CHANGELOG`, `docs/src/*.md`) → `MIT OR CC-BY-4.0`; trivial config/tooling (`Makefile`, `.gitignore`, CI YAML, `Project.toml`, …) → `MIT OR CC0-1.0`. Any license used must have its text in `LICENSES/`.

[JuliaQuantumControl]: https://github.com/JuliaQuantumControl
