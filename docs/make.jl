using QuantumPropagators
using Documenter

DocMeta.setdocmeta!(QuantumPropagators, :DocTestSetup, :(using QuantumPropagators); recursive=true)

makedocs(;
    modules=[QuantumPropagators],
    authors="Michael Goerz <mail@michaelgoerz.net> and contributors",
    repo="https://github.com/quantumcontrol-jl/QuantumPropagators.jl/blob/{commit}{path}#{line}",
    sitename="QuantumPropagators.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://quantumcontrol-jl.github.io/QuantumPropagators.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/quantumcontrol-jl/QuantumPropagators.jl",
)
