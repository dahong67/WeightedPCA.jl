using WeightedPCA
using Documenter

DocMeta.setdocmeta!(WeightedPCA, :DocTestSetup, :(using WeightedPCA); recursive=true)

makedocs(;
    modules=[WeightedPCA],
    authors="David Hong <dahong67@wharton.upenn.edu> and contributors",
    repo="https://github.com/dahong67/WeightedPCA.jl/blob/{commit}{path}#{line}",
    sitename="WeightedPCA.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://dahong67.github.io/WeightedPCA.jl",
        edit_link="master",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/dahong67/WeightedPCA.jl",
    devbranch="master",
)
