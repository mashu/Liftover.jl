using Documenter
using Documenter.Remotes: GitHub
using Liftover

makedocs(
    modules = [Liftover],
    sitename = "Liftover.jl",
    authors = "Mateusz Kaduk",
    repo = GitHub("mashu", "Liftover.jl"),
    format = Documenter.HTML(
        prettyurls = get(ENV, "CI", nothing) == "true",
        size_threshold_ignore = ["api.md"],
        edit_link = "main",
    ),
    pages = [
        "Home" => "index.md",
        "Getting started" => "getting_started.md",
        "Coordinates" => "coordinates.md",
        "Examples" => "examples.md",
        "API" => "api.md",
    ],
    checkdocs = :exports,
)

deploydocs(
    repo = "github.com/mashu/Liftover.jl.git",
    devbranch = "main",
    push_preview = true,
)
