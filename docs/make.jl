using RadarData
using Documenter

DocMeta.setdocmeta!(RadarData, :DocTestSetup, :(using RadarData); recursive=true)

makedocs(;
    modules=[RadarData],
    authors="Sam Pine <sam.pine@matrixresearch.com> and contributors",
    sitename="RadarData.jl",
    format=Documenter.HTML(;
        canonical="https://Programs/alg/julia.gitlab.io/RadarData.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)
