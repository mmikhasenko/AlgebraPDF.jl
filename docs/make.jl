using Pkg
cd(joinpath(@__DIR__, ".."))
Pkg.activate("docs")
Pkg.develop(path=".")

using Documenter
using AlgebraPDF

function fix_literate_output(content)
    content = replace(content, "EditURL = \"@__REPO_ROOT_URL__/\"" => "")
    return content
end


gen_content_dir = joinpath(@__DIR__, "src")
tutorial_src = joinpath(@__DIR__, "src", "tutorial_lit.jl")
Literate.markdown(tutorial_src, gen_content_dir, name = "tutorial", documenter = true, credit = true, postprocess = fix_literate_output)
Literate.notebook(tutorial_src, gen_content_dir, execute = false, name = "algebrapdf_tutorial", documenter = true, credit = true)
Literate.script(tutorial_src, gen_content_dir, keep_comments = false, name = "algebrapdf_tutorial", documenter = true, credit = false)


makedocs(
    sitename="AlgebraPDF",
    modules=[AlgebraPDF],
    format = Documenter.HTML(
        prettyurls = false,
        canonical = "https://mmikhasenko.github.io/algebrapdf.jl/stable/"
    ),
    authors = "Mikhail Mikhasenko",
    pages=[
        "Home" => "index.md",
        "Tutorial" => "tutorial.md",
        "API" => "api.md"
    ],
    doctest = true,
    linkcheck = false
)

deploydocs(
    repo="github.com/mmikhasenko/AlgebraPDF.jl.git",
    target="site",
)
