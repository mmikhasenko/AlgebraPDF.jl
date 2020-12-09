# AlgebraPDF

[![Build Status](https://travis-ci.com/mmikhasenko/AlgebraPDF.jl.svg?branch=master)](https://travis-ci.com/mmikhasenko/AlgebraPDF.jl)
[![Build Status](https://ci.appveyor.com/api/projects/status/github/mmikhasenko/AlgebraPDF.jl?svg=true)](https://ci.appveyor.com/project/mmikhasenko/AlgebraPDF-jl)
[![Codecov](https://codecov.io/gh/mmikhasenko/AlgebraPDF.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/mmikhasenko/AlgebraPDF.jl)

Basic functionality:
```julia
using AlgebraPDF
using Plots
pyplot()
theme(:wong)

# create a model
snl = pdf(@. (x;p) -> exp(-(x-p.μ)^2/(2*p.σ^2)); p = (μ=1.4, σ=0.15), lims=(0, 3))
bkg = pdf(@. (x;p) -> sqrt(x)*exp(-p.α*x); p = (α=1.3,), lims=(0, 3))
bkg *= (fb=2.5,)
pdf_sum = snl + bkg

# generating
const data = generate(1000, pdf_sum);

# fitting
pfr = v2p(fit_llh(data, pdf_sum; init_pars=p2v(pdf_sum)), pdf_sum)

# plotting
let
  plot()
  plot!(x->pdf_sum(x; p=pfr), lims(pdf_sum)..., lab="fit")
  plot!(x->snl(x; p=pfr, norm_according_to=pdf_sum), lims(pdf_sum)..., lab="signal")
  plot!(x->bkg(x; p=pfr, norm_according_to=pdf_sum), lims(pdf_sum)..., lab="background")
  stephist!(data, norm=true, c=:black, bins=50, lab="data")
end
```
![example](plots/gaus.background.png)


# Todo:
 - fixing parameters
 - convolution
 - perhaps, add basic function: `gaus`, `BW`, `polN`
 - add documentation
