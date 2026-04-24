push!(LOAD_PATH, dirname(@__DIR__))

using Documenter
using RheologyCalculator

const CI = get(ENV, "CI", nothing) == "true"

@info "Making documentation..."
makedocs(;
    sitename = "RheologyCalculator.jl",
    authors  = "Albert de Montserrat and Boris Kaus",
    modules  = [RheologyCalculator],
    format   = Documenter.HTML(;
        prettyurls = CI,
        disable_git = !CI,
        edit_link = nothing,
    ),
    warnonly = Documenter.except(:footnote),
    checkdocs = :exports,
    pages    = [
        "Home" => "index.md",
        "Composites" => "composites.md",
        "Rheology" => "rheology.md",
        "Solving" => "stress.md",
        "API" => "api.md",
    ],
)

deploydocs(; repo="https://github.com/albert-de-montserrat/RheologyCalculator.jl", devbranch="main")
