using PkgTemplates

t = Template(;
           user="mmikhasenko",
           host="github.com",
           license="MIT",
           authors="Misha Mikhasenko",
           julia_version=v"1.3",
           dir=joinpath(DEPOT_PATH[1], "dev"),
           ssh=true,
           plugins=[
               TravisCI(),
               Codecov(),
               AppVeyor(),
           ],
       )

generate("AlgebraPDF", t)
