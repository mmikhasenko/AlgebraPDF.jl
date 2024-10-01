# AlgebraPDF

[![Stable Documentation](https://img.shields.io/badge/docs-stable-blue.svg)](https://mmikhasenko.github.io/AlgebraPDF.jl/stable)
[![In development documentation](https://img.shields.io/badge/docs-dev-blue.svg)](https://mmikhasenko.github.io/AlgebraPDF.jl/dev)
[![Build Status](https://github.com/mmikhasenko/AlgebraPDF.jl/workflows/Test/badge.svg)](https://github.com/mmikhasenko/AlgebraPDF.jl/actions)
[![Test workflow status](https://github.com/mmikhasenko/AlgebraPDF.jl/actions/workflows/Test.yml/badge.svg?branch=main)](https://github.com/mmikhasenko/AlgebraPDF.jl/actions/workflows/Test.yml?query=branch%3Amain)
[![Lint workflow Status](https://github.com/mmikhasenko/AlgebraPDF.jl/actions/workflows/Lint.yml/badge.svg?branch=main)](https://github.com/mmikhasenko/AlgebraPDF.jl/actions/workflows/Lint.yml?query=branch%3Amain)
[![Docs workflow Status](https://github.com/mmikhasenko/AlgebraPDF.jl/actions/workflows/Docs.yml/badge.svg?branch=main)](https://github.com/mmikhasenko/AlgebraPDF.jl/actions/workflows/Docs.yml?query=branch%3Amain)

[![Coverage](https://codecov.io/gh/mmikhasenko/AlgebraPDF.jl/branch/main/graph/badge.svg)](https://codecov.io/gh/mmikhasenko/AlgebraPDF.jl)
[![DOI](https://zenodo.org/badge/DOI/FIXME)](https://doi.org/FIXME)

## Basic functionality

* Attach default values of parameters to a function
* Update, fix, release parameters
* constructing a complex model object from set of function, e.g. `f₁ + f₂`, `abs2(f)`, or `log(f)`.
* On-fly normalization
* construction of mixed models in the form `f₁ PDF₁ + f₂ PDF₂ + f₃ PDF₃`.
* construction of likelihood function and extended likelihood function
* plotting recipes

Current implementation is limited to immutable objects.

## How to Cite

If you use AlgebraPDF.jl in your work, please cite using the reference given in [CITATION.cff](https://github.com/mmikhasenko/AlgebraPDF.jl/blob/main/CITATION.cff).
