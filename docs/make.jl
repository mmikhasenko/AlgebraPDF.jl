using Pkg
cd(joinpath(@__DIR__, ".."))
Pkg.activate("docs")  # Activate the documentation environment
Pkg.develop(path=".")  # Make sure the package is available in the environment

using Documenter
using AlgebraPDF  # Now you can use your package

makedocs(
    sitename="AlgebraPDF Documentation",
    modules=[AlgebraPDF],
    pages=[
        "Home" => "index.md",
        "Tutorial" => "tutorial.md"
    ]
)

deploydocs(
    repo="github.com/mmikhasenko/AlgebraPDF.jl.git",
    target="site",
)
