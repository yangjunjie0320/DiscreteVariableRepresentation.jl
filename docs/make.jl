using DiscreteVariableRepresentation
using Documenter

DocMeta.setdocmeta!(DiscreteVariableRepresentation, :DocTestSetup, :(using DiscreteVariableRepresentation); recursive=true)

makedocs(;
    modules=[DiscreteVariableRepresentation],
    authors="yangjunjie",
    sitename="DiscreteVariableRepresentation.jl",
    format=Documenter.HTML(;
        canonical="https://yangjunjie.github.io/DiscreteVariableRepresentation.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/yangjunjie/DiscreteVariableRepresentation.jl",
    devbranch="main",
)
