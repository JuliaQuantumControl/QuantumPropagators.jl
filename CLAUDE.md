# CLAUDE.md

This file provides guidance to Claude Code (claude.ai/code) when working with code in this repository.

## Development Commands

### Testing
- `make test` - Run the full test suite
- `julia --project=test -e 'include("test/runtests.jl")'` - Run tests directly

### Development Environment
- `make devrepl` - Start interactive REPL with development environment (recommended for development)
- `julia -i --banner=no devrepl.jl` - Alternative way to start development REPL

### Documentation
- `make docs` - Build documentation

### Code Quality
- `make codestyle` - Apply JuliaFormatter to entire project
- `julia --project=test -e 'using JuliaFormatter; format(".", verbose=true)'` - Format code directly

### Cleanup
- `make clean` - Clean build/doc/testing artifacts
- `make distclean` - Restore to clean checkout state

## Project Architecture

This is a Julia package for quantum system time propagation within the JuliaQuantumControl ecosystem.

### Core Structure
- **src/QuantumPropagators.jl** - Main module definition with submodule includes
- **src/interfaces/** - Interface definitions for operators, states, generators, propagators, etc.
- **Submodules** organized by functionality:
  - `Arnoldi` (arnoldi.jl) - Arnoldi iteration methods
  - `SpectralRange` (specrad.jl) - Spectral radius calculations
  - `Cheby` (cheby.jl) - Chebychev propagation methods
  - `Newton` (newton.jl) - Newton propagation methods
  - `ExpProp` (expprop.jl) - Matrix exponential propagation
  - `Storage` (storage.jl) - Memory management utilities
  - `Shapes` (shapes.jl) - Pulse shape functions
  - `Controls` (controls.jl) - Control parameter handling
  - `Amplitudes` (amplitudes.jl) - Time-dependent amplitude functions
  - `Generators` (generators.jl) - Hamiltonian/Liouvillian generators

### Key Propagation Methods
- **Chebyshev propagation** (cheby_propagator.jl) - Polynomial expansion method
- **Newton propagation** (newton_propagator.jl) - Newton interpolation method
- **Matrix exponential** (exp_propagator.jl) - Direct exponentiation
- **ODE integration** (ode_function.jl) - Interface to OrdinaryDiffEq.jl via extensions

### High-Level Interface
- `propagate()` function for single time evolution
- `propagate_sequence()` for sequential propagations with different parameters

### Extensions
The package uses Julia's extension system for optional dependencies:
- `QuantumPropagatorsODEExt` - OrdinaryDiffEq.jl integration
- `QuantumPropagatorsRecursiveArrayToolsExt` - RecursiveArrayTools.jl support
- `QuantumPropagatorsStaticArraysExt` - StaticArrays.jl optimization

### Testing Framework
Uses SafeTestsets for isolated test execution. Tests are comprehensive and include:
- Interface validation tests
- Individual propagator method tests
- Integration tests for complete propagation workflows

### Development Environment
The project uses a sophisticated development setup:
- Development REPL (devrepl.jl) with Revise.jl for hot reloading
- Automatic dependency management via installorg.jl script
- Integrated documentation building and serving
- Code formatting with JuliaFormatter

### Documentation
Uses Documenter.jl with comprehensive API documentation and examples. Documentation includes detailed method explanations.

## Docstrings

Each Julia function that is not explicitly private or has a name starting with an underscore must have a docstring. The format for the docstring is to have a "title" / description of the function in the first line, less than 80 characters, followed by a fenced code block (```` ```julia...``` ````; do not use an indented code block). The code block should be valid Julia code that illustrates a default call to the function. If the function returns something, the code block should be an assignment with instructive names for the returned variables. The code block should be readable as the start of an active-voice sentence that continues after the code block. Following the paragraph started by the code block, use section `#Arguments`, `# Keyword Arguments`, `# Returns`, etc, as necessary. Be as concise as possible. For example, if the code block paragraph already fully explains the return type, an explicit `# Returns` section is unnecessary. Use Unicode as much as possible for any math in the docstring.

## General Guidelines

* Make sure to only use explicit imports in Julia code, and that there are no imported functions or constants that are not actually used.

* When adding a new dependency to any `Project.toml` file, run `make distclean`, and then `make test/Manifest.toml`, `make docs/Manifest.toml`, etc. to recreate manifest files as necessary.

* Never commit any changes or ask to commit. I will always create git commits manually.
