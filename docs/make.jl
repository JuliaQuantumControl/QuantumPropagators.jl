using Pkg

using Documenter
using QuantumPropagators
using QuantumPropagators: AbstractPropagator, set_t!, set_state!
import OrdinaryDiffEq  # ensure ODE extension is loaded
using DocumenterCitations
using DocumenterInterLinks


PROJECT_TOML = Pkg.TOML.parsefile(joinpath(@__DIR__, "..", "Project.toml"))
VERSION = PROJECT_TOML["version"]
NAME = PROJECT_TOML["name"]
AUTHORS = join(PROJECT_TOML["authors"], ", ") * " and contributors"
GITHUB = "https://github.com/JuliaQuantumControl/QuantumPropagators.jl"

DEV_OR_STABLE = "stable/"
if endswith(VERSION, "dev")
    DEV_OR_STABLE = "dev/"
end


function org_inv(pkgname)
    objects_inv = joinpath(
        @__DIR__, "..", "..", "$pkgname.jl", "docs", "build", "objects.inv"
    )
    if isfile(objects_inv)
        return (
            "https://juliaquantumcontrol.github.io/$pkgname.jl/dev/",
            objects_inv,
        )
    else
        return "https://juliaquantumcontrol.github.io/$pkgname.jl/$DEV_OR_STABLE"
    end
end


links = InterLinks(
    "Julia" => (
        "https://docs.julialang.org/en/v1/",
        "https://docs.julialang.org/en/v1/objects.inv",
        joinpath(@__DIR__, "src", "inventories", "Julia.toml"),
    ),
    "TimerOutputs" => (
        "https://github.com/KristofferC/TimerOutputs.jl",
        joinpath(@__DIR__, "src", "inventories", "TimerOutputs.toml"),
    ),
    "QuantumControl" => org_inv("QuantumControl"),
    "StaticArrays" => "https://juliaarrays.github.io/StaticArrays.jl/stable/",
    "ComponentArrays" => "https://jonniedie.github.io/ComponentArrays.jl/stable/",
    "RecursiveArrayTools" => "https://docs.sciml.ai/RecursiveArrayTools/stable/",
    "qutip" => "https://qutip.readthedocs.io/en/qutip-5.0.x/",
)

externals = ExternalFallbacks()

println("Starting makedocs")

include("generate_api.jl")

bib = CitationBibliography(joinpath(@__DIR__, "src", "refs.bib"); style=:numeric)

warnonly = [:linkcheck,]
if get(ENV, "DOCUMENTER_WARN_ONLY", "0") == "1"  # cf. test/init.jl
    warnonly = true
end


makedocs(;
    plugins=[bib, links, externals],
    authors=AUTHORS,
    sitename="QuantumPropagators.jl",
    # Link checking is disabled in REPL, see `devrepl.jl`.
    linkcheck=(get(ENV, "DOCUMENTER_CHECK_LINKS", "1") != "0"),
    warnonly,
    doctest=false,  # doctests run as part of test suite
    format=Documenter.HTML(;
        prettyurls=true,
        canonical="https://juliaquantumcontrol.github.io/QuantumPropagators.jl",
        assets=[
            "assets/citations.css",
            asset(
                "https://juliaquantumcontrol.github.io/QuantumControl.jl/dev/assets/topbar/topbar.css"
            ),
            asset(
                "https://juliaquantumcontrol.github.io/QuantumControl.jl/dev/assets/topbar/topbar.js"
            ),
        ],
        footer="[$NAME.jl]($GITHUB) v$VERSION docs powered by [Documenter.jl](https://github.com/JuliaDocs/Documenter.jl).",
        size_threshold_ignore=["api/quantumpropagators.md",]
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
    ],
)

println("Finished makedocs")

deploydocs(; repo="github.com/JuliaQuantumControl/QuantumPropagators.jl")
