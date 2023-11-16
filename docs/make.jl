using Pkg
cd(joinpath(@__DIR__, ".."))
Pkg.activate("docs")
Pkg.develop(PackageSpec(path=pwd()))  # adds the package that is being documented
Pkg.instantiate()  # download & install dependences specific for the documentation

using Documenter
using AlgebraPDF
using Literate

function fix_literate_output(content)
    content = replace(content, "EditURL = \"@__REPO_ROOT_URL__/\"" => "")
    return content
end


gen_content_dir = joinpath(@__DIR__, "src")
tutorial_src = joinpath(@__DIR__, "src", "tutorial_lit.jl")
Literate.markdown(tutorial_src, gen_content_dir, name="tutorial", documenter=true, credit=true, postprocess=fix_literate_output)
Literate.notebook(tutorial_src, gen_content_dir, execute=false, name="algebrapdf_tutorial", documenter=true, credit=true)
Literate.script(tutorial_src, gen_content_dir, keep_comments=false, name="algebrapdf_tutorial", documenter=true, credit=false)


makedocs(
    sitename="AlgebraPDF",
    modules=[AlgebraPDF],
    format=Documenter.HTML(
        prettyurls=true,
        canonical="https://mmikhasenko.github.io/AlgebraPDF.jl/stable/"
    ),
    authors="Mikhail Mikhasenko",
    pages=[
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "API" => "api.md"
    ],
    doctest=false,
    linkcheck=false
)

deploydocs(
    repo="github.com/mmikhasenko/AlgebraPDF.jl.git",
    push_preview=true,
)
