using RadarData
using Documenter

DocMeta.setdocmeta!(RadarData, :DocTestSetup, :(using RadarData); recursive = true)

makedocs(;
    modules = [RadarData],
    authors = "Sam Pine <sam.pine@matrixresearch.com> and contributors",
    sitename = "RadarData.jl",
    format = Documenter.HTML(;
        canonical = "https://git.matrixresearch.com/Programs/alg/julia/radardata",
        edit_link = "main",
        assets = String[]
    ),
    pages = [
        "Home" => "index.md"
    ]
)
