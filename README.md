# AlgebraPDF

[![Build Status](https://github.com/mmikhasenko/AlgebraPDF.jl/workflows/CI/badge.svg)](https://github.com/mmikhasenko/AlgebraPDF.jl/actions)
[![Codecov](https://codecov.io/gh/mmikhasenko/AlgebraPDF.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/mmikhasenko/AlgebraPDF.jl)

Basic functionality:
 * constructing a pdf object from arbitraty function:
    - required function format `myfunc(x;p)` where `p` is a NamedTuple of parameters
    - the normalization is computed automatically using `QuadGK.jl` in the giben `lims`
    - algebra of functions with parameters: `*`, `+`, `-`, `/`
 * construction of mixed models in the form `f₁ PDF₁ + f₂ PDF₂ + (1-f₁-f₂) PDF₃`.
 * fitting data distribution using the maximum likelihood (`Optim.jl`)
 * plotting recipies

```julia
using AlgebraPDF
using Plots
theme(:wong)

# create a model
snl = aGauss((μ=1.4, σ=0.15), (0, 3))
bkg = pdf(@. (x;p) -> sqrt(x)*exp(-p.α*x); p = (α=1.3,), lims=(0, 3))
bkg_f = bkg * (fb=2.5,)
pdf_sum = snl + bkg_f

# generating
const data = generate(1000, pdf_sum);

# fitting
fr = fit_llh(data, pdf_sum; init_pars=p2v(pdf_sum))
pfr = v2p(minimizer(fr), pdf_sum) # NamedTuple (parameter = value, ...)
efr = v2p(errors(fr), pdf_sum) # NamedTuple (parameter = error, ...)
mfr = v2p(measurements(fr), pdf_sum) # NamedTuple (parameter = value ± error, ...)

# plotting
let
  plot()
  plot!(x->pdf_sum(x; p=pfr), lims(pdf_sum)..., lab="fit")
  plot!(x->snl(x; p=pfr, norm_according_to=pdf_sum), lims(pdf_sum)..., lab="signal")
  plot!(x->bkg(x; p=pfr, norm_according_to=pdf_sum), lims(pdf_sum)..., lab="background")
  stephist!(data, norm=true, c=:black, bins=50, lab="data")
end

# Alternatively one can use a mixed model
model = MixedModel([snl, bkg], (fS = 0.5,))
fr = fit_llh(data, model)
pfr = v2p(minimizer(fr), pdf_sum)
fixed_model = fixpars(model, pfr)
# 
let
  Nd = length(data)
  bins = range(lims(fixed_model)..., length = 40)
  Ns = scaletobinneddata(Nd, bins)
  #
  plot(pdf_sum, Ns, lab="fit")
  plot!(pdf_sum.components[1], fractionvalues(fixed_model)[1]*Ns, lab="signal")
  plot!(pdf_sum.components[2], fractionvalues(fixed_model)[2]*Ns, lab="background")
  stephist!(data, c=:black, bins=bins, lab="data")
end
```
![example](plots/gaus.background.png)





```julia
# implementation with NAMES of parameters build into the funciton call
struct BW1{P} <: AbstractFunctionWithParameters
    p::P
end
import AlgebraPDF:func
func(bw::BW1, x::Number; p=pars(bw)) = p.m*p.Γ/(p.m^2-x^2-1im*p.m*p.Γ)

bw = BW1((m=0.77, Γ=0.15))
# 
bw(0.77)  # 0.0 + 0.9999999999999999im
bw(0.77; p=(m=0.8, Γ=0.15))  # 0.34010473926206014 + 0.866508889839641im
```

Slightly better implementation where only the order of the arguments are fixed, while the names are determined when the instance is created.
```julia
# implementation with ORDER of parameters build into the funciton call
struct BW2{P} <: AbstractFunctionWithParameters
    p::P
end
import AlgebraPDF:func
function func(bw::BW2, x::Number; p=pars(bw))
    m,Γ = (getproperty(p,s) for s in keys(bw.p))
    m*Γ/(m^2-x^2-1im*m*Γ)
end

# same function with different names
bw_i = BW2((m_i=1.1, Γ_i=0.2))
bw_j = BW2((m_j=1.1, Γ_j=0.2))
bw_k = BW2((m_k=1.1, Γ_k=0.2))
```