using QuantumPropagators
using Documenter

DocMeta.setdocmeta!(QuantumPropagators, :DocTestSetup, :(using QuantumPropagators); recursive=true)

makedocs(;
    modules=[QuantumPropagators],
    authors="Michael Goerz <mail@michaelgoerz.net> and contributors",
    repo="https://github.com/JuliaQuantumControl/QuantumPropagators.jl/blob/{commit}{path}#{line}",
    sitename="QuantumPropagators.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://juliaquantumcontrol.github.io/QuantumPropagators.jl",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "Overview" => "overview.md",
        "Background" => "background.md",
        hide("Examples" => "examples/index.md", [
            # "Example 1" => joinpath("examples", "1.md"),
            # "Example 2" => joinpath("examples", "2.md"),
           ]
        ),
        "Howtos" => "howto.md",
        "Benchmarks" => "benchmarks.md",
        "API" => "api.md",
    ],
)

deploydocs(;
    repo="github.com/JuliaQuantumControl/QuantumPropagators.jl",
)
