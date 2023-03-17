using QuantumPropagators
using Documenter
using Pkg
using QuantumCitations

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

include("generate_api.jl")

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style=:numeric)

makedocs(
    bib;
    authors=AUTHORS,
    sitename="QuantumPropagators.jl",
    format=Documenter.HTML(;
        prettyurls=true,
        canonical="https://juliaquantumcontrol.github.io/QuantumPropagators.jl",
        assets=String["assets/citations.css"],
        footer="[$NAME.jl]($GITHUB) v$VERSION docs powered by [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl)."
    ),
    pages=[
        "Home" => "index.md",
        "Overview" => "overview.md",
        "Dynamical Generators" => "generators.md",
        "Propagation Methods" => "methods.md",
        "Storage" => "storage.md",
        hide("Examples" => "examples/index.md", [
        # "Example 1" => joinpath("examples", "1.md"),
        # "Example 2" => joinpath("examples", "2.md"),
        ]),
        "Howtos" => "howto.md",
        "Benchmarks" => "benchmarks.md",
        "API" => "api/quantumpropagators.md",
        "References" => "references.md",
    ]
)

println("Finished makedocs")

deploydocs(; repo="github.com/JuliaQuantumControl/QuantumPropagators.jl")
