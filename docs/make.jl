using Documenter, UKrig

makedocs(;
    modules=[UKrig],
    format=Documenter.HTML(),
    pages=[
        "Home" => "index.md",
    ],
    repo="https://github.com/hango1996/UKrig.jl/blob/{commit}{path}#L{line}",
    sitename="UKrig.jl",
    authors="Han Chen",
    assets=String[],
)

deploydocs(;
    repo="github.com/hango1996/UKrig.jl",
)
