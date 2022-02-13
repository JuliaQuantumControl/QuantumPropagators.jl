# init file for "make devrepl"
using Revise
using Plots
unicodeplots()
using JuliaFormatter
using QuantumControlBase.TestUtils: test
using LiveServer: LiveServer, serve, servedocs
include(joinpath(@__DIR__, "clean.jl"))


println("""
*******************************************************************************
DEVELOPMENT REPL

Revise, JuliaFormatter, LiveServer, Plots with unicode backend are active.

* `include("test/runtests.jl")` – Run the entire test suite
* `test()` – Run the entire test suite in a subprocess with coverage
* `test(genhtml=true)` – Generate an HTML coverage report
* `Pkg.test("QuantumControlBase", coverage=true)` –
  Run upstream QuantumControlBase tests for additional coverage
* `include("docs/make.jl")` – Generate the documentation
* `format(".")` – Apply code formatting to all files
* `servedocs([port=8000, verbose=false])` –
  Build and serve the documentation. Automatically recompile and redisplay on
  changes
* `clean()` – Clean up build/doc/testing artifacts
* `distclean()` – Restore to a clean checkout state
*******************************************************************************
""")
