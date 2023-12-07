using QuantumPropagators
using QuantumPropagators: AbstractPropagator, set_t!, set_state!
using Documenter
using Pkg
using OrdinaryDiffEq  # ensure ODE extension is loaded
using DocumenterCitations

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

warnonly = [:linkcheck,]
if get(ENV, "DOCUMENTER_WARN_ONLY", "0") == "1"  # cf. test/init.jl
    warnonly = true
end


makedocs(;
    plugins=[bib],
    authors=AUTHORS,
    sitename="QuantumPropagators.jl",
    # Link checking is disabled in REPL, see `devrepl.jl`.
    linkcheck=(get(ENV, "DOCUMENTER_CHECK_LINKS", "1") != "0"),
    warnonly,
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
        "Expectation Values" => "storage.md",
        "Howtos" => "howto.md",
        hide("Benchmarks" => "benchmarks.md", [joinpath("benchmarks", "profiling.md")]),
        "API" => "api/quantumpropagators.md",
        "References" => "references.md",
    ]
)

println("Finished makedocs")

deploydocs(; repo="github.com/JuliaQuantumControl/QuantumPropagators.jl")
