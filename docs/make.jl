using Documenter
using PDESystemLibrary

makedocs(;
    modules = [PDESystemLibrary],
    checkdocs = :exports,
    sitename = "PDESystemLibrary.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://docs.sciml.ai/PDESystemLibrary/stable/",
        edit_link = "master",
    ),
    pages = [
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo = "github.com/SciML/PDESystemLibrary.jl",
    devbranch = "master",
)
