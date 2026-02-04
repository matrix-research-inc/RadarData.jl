using RadarData
using Documenter

DocMeta.setdocmeta!(RadarData, :DocTestSetup, :(using RadarData); recursive = true)

makedocs(;
    modules = [RadarData],
    authors = "Sam Pine <sam.pine@matrixresearch.com> and contributors",
    sitename = "RadarData.jl",
    remotes = nothing,
    format = Documenter.HTML(;
        canonical = "https://git.matrixresearch.com/Programs/alg/julia/radardata",
        edit_link = "main",
        prettyurls = get(ENV, "CI", "false") == "true",
        assets = String[]
    ),
    warnonly = [:missing_docs, :cross_references, :docs_block],
    pages = [
        "Home" => "index.md"
    ]
)
