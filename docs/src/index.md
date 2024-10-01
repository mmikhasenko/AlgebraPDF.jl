```@meta
CurrentModule = AlgebraPDF
```

# AlgebraPDF

Documentation for [AlgebraPDF](https://github.com/mmikhasenko/AlgebraPDF.jl). AlgebraPDF.jl is a tool to construct and operate with complex parametric probability density functions (PDFs).

Basic functionality:

* Attach default values of parameters to a function
* Update, fix, release parameters
* constructing a complex model object from set of function:
    * algebra of functions with parameters, e.g. `f₁ + f₂`, `abs2(f)`, or `log(f)`.
* On-fly normalization
* construction of mixed models in the form `f₁ PDF₁ + f₂ PDF₂ + f₃ PDF₃`.
* construction of likelihood function and extended likelihood function
* plotting recipes

!!! note

    `AlgebraPDF.jl` is experimental, feel free to share experience working with the package. Submit an issue with suggestions.

## Table of contents

```@contents
Pages = [
    "tutorial.md",
    "api.md",
    "developing.md",
    "license.md",
]
Depth = 1
```

## Citing AlgebraPDF.jl

When using `AlgebraPDF.jl` for research, teaching or similar,
please drop a link for your research to the issues.
It would give me encouragement to create a citable link on Zenodo, write a JOSS paper. I would also be happy to mention your research in the Documentation website.

## Learning Julia

The [Julia website](https://julialang.org/) provides many [links to introductory videos and written tutorials](https://julialang.org/learning/), e.g. ["Intro to Julia"](https://www.youtube.com/watch?v=fMa1qSg_LxA),
[Think Julia: How to Think Like a Computer Scientist](https://benlauwens.github.io/ThinkJulia.jl/latest/book.html)
and ["The Fast Track to Julia"](https://juliadocs.github.io/Julia-Cheat-Sheet/). If you are familiar with MATLAB or Python, you may also want to take a look at the ["MATLAB–Python–Julia cheatsheet"](https://cheatsheets.quantecon.org/).

The in-depth article [Why Numba and Cython are not substitutes for Julia](http://www.stochasticlifestyle.com/why-numba-and-cython-are-not-substitutes-for-julia/) explains how Julia addresses several fundamental challenges inherent to scientific high-performance computing.

## Installation

```julia
using Pkg
Pkg.add("AlgebraPDF")
```

## Usage

Provide examples of how to use your package.
