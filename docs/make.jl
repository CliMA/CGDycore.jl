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
        "Home" => "index.md",
        "Mass Matrix" => "MassMatrix.md",
        "Stiffness Matrix" => "StiffMatrix.md",
        "Examples" => "Examples.md",
    ],
)

deploydocs(; repo = "github.com/CliMA/CGDycore.jl", push_preview = true)