# init file for "make devrepl"
using Pkg
Pkg.develop(PackageSpec(path=pwd()))
using Revise
println("""
*******************************************************************************
DEVELOPMENT REPL

Revise is active

Run

    include("test/runtests.jl")

for running the entire test suite. Alternatively, run e.g.

    include("test/test_propagate.jl")

for running an individual test file. Run

    include("docs/make.jl")

to generate the documentation
*******************************************************************************
""")
