using ChebParticleMesh
using Documenter

DocMeta.setdocmeta!(ChebParticleMesh, :DocTestSetup, :(using ChebParticleMesh); recursive=true)

makedocs(;
    modules=[ChebParticleMesh],
    authors="xzgao <xz.gao@connect.ust.hk> and contributors",
    repo="https://github.com/Xuanzhao Gao/ChebParticleMesh.jl/blob/{commit}{path}#{line}",
    sitename="ChebParticleMesh.jl",
    format=Documenter.HTML(;
        prettyurls=get(ENV, "CI", "false") == "true",
        canonical="https://Xuanzhao Gao.github.io/ChebParticleMesh.jl",
        edit_link="main",
        assets=String[],
    ),
    pages=[
        "Home" => "index.md",
    ],
)

deploydocs(;
    repo="github.com/Xuanzhao Gao/ChebParticleMesh.jl",
    devbranch="main",
)
