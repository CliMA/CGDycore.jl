using CGDycore
using Documenter

makedocs(;
    modules = [
        CGDycore,
    ],
    authors = "Oswald Knoth",
    sitename = "CGDycore.jl",
    format = Documenter.HTML(;
        prettyurls = !isempty(get(ENV, "CI", "")),
        collapselevel = 1,
    ),
    checkdocs = :exports,
    pages = [
        "Mass Matrix" => "MassMatrix.md",

    ],
)

deploydocs(; repo = "github.com/CliMA/CGDycore.jl", push_preview = true)
