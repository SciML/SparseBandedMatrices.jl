using SparseBandedMatrices
using Documenter

cp("./docs/Manifest.toml", "./docs/src/assets/Manifest.toml", force = true)
cp("./docs/Project.toml", "./docs/src/assets/Project.toml", force = true)

makedocs(;
    modules = [SparseBandedMatrices],
    authors = "Chris Rackauckas <accounts@chrisrackauckas.com> and contributors",
    repo = "https://github.com/SciML/SparseBandedMatrices.jl/blob/{commit}{path}#{line}",
    sitename = "SparseBandedMatrices.jl",
    format = Documenter.HTML(;
        prettyurls = get(ENV, "CI", "false") == "true",
        canonical = "https://docs.sciml.ai/SparseBandedMatrices/stable/",
        edit_link = "main",
        assets = String[],
    ),
    pages = [
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo = "github.com/SciML/SparseBandedMatrices.jl",
    devbranch = "main",
)
