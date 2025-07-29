using Documenter, RheologyCalculator
push!(LOAD_PATH, "../src/")

@info "Making documentation..."
makedocs(;
    sitename = "RheologyCalculator.jl",
    authors  = "Albert de Montserrat and Boris Kaus",
    modules  = [RheologyCalculator],
    format   = Documenter.HTML(; prettyurls=get(ENV, "CI", nothing) == "true"), # easier local build
    warnonly = Documenter.except(:footnote),
    pages    = [
        "Home" => "index.md",
        "Rheology" => "rheology.md",
    ],
)

deploydocs(; repo="https://github.com/albert-de-montserrat/RheologyCalculator.jl", devbranch="main")
