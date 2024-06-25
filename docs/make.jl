using RGRDMT
using Documenter

DocMeta.setdocmeta!(RGRDMT, :DocTestSetup, :(using RGRDMT); recursive=true)

makedocs(;
    modules=[RGRDMT],
    authors="Yusheng Zhao <yushengzhao2020@outlook.com> and contributors",
    sitename="RGRDMT.jl",
    format=Documenter.HTML(;
        canonical="https://exAClior.github.io/RGRDMT.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
        "API" => "api.md"
    ],
)

deploydocs(;
    repo="github.com/exAClior/RGRDMT.jl",
)

