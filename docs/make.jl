using QuantumPropagators
using Documenter
using Pkg

DocMeta.setdocmeta!(
    QuantumPropagators,
    :DocTestSetup,
    :(using QuantumPropagators);
    recursive=true
)

PROJECT_TOML = Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))
VERSION = PROJECT_TOML["version"]
NAME = PROJECT_TOML["name"]
AUTHORS = join(PROJECT_TOML["authors"], ", ") * " and contributors"
GITHUB = "https://github.com/JuliaQuantumControl/QuantumPropagators.jl"

println("Starting makedocs")

makedocs(;
    authors=AUTHORS,
    sitename="QuantumPropagators.jl",
    modules=[QuantumPropagators],
    repo="https://github.com/JuliaQuantumControl/QuantumPropagators.jl/blob/{commit}{path}#{line}",
    format=Documenter.HTML(;
        prettyurls=true,
        canonical="https://juliaquantumcontrol.github.io/QuantumPropagators.jl",
        assets=String[],
        footer="[$NAME.jl]($GITHUB) v$VERSION docs powered by [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl)."
    ),
    pages=[
        "Home" => "index.md",
        "Overview" => "overview.md",
        "Background" => "background.md",
        hide("Examples" => "examples/index.md", [
        # "Example 1" => joinpath("examples", "1.md"),
        # "Example 2" => joinpath("examples", "2.md"),
        ]),
        "Howtos" => "howto.md",
        "Benchmarks" => "benchmarks.md",
        "API" => "api.md",
        "History" => "history.md",
    ]
)

println("Finished makedocs")

deploydocs(; repo="github.com/JuliaQuantumControl/QuantumPropagators.jl")
