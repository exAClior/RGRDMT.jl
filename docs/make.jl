using RGRDMT
using Documenter
using DocumenterCitations


DocMeta.setdocmeta!(RGRDMT, :DocTestSetup, :(using RGRDMT); recursive=true)

bib = CitationBibliography(joinpath(@__DIR__,"src/reference.bib"),style=:authoryear)

makedocs(;
    plugins = [bib],
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
        "Theory" => "theory.md",
        "Implementation" => "implementation.md",
        "API" => "api.md",
        "Suggested Readings and References" => "reference.md"
    ],
)

deploydocs(;
    repo="github.com/exAClior/RGRDMT.jl",
)

