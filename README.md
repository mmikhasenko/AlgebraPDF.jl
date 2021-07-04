# AlgebraPDF

[![Build Status](https://github.com/mmikhasenko/AlgebraPDF.jl/workflows/CI/badge.svg)](https://github.com/mmikhasenko/AlgebraPDF.jl/actions)
[![Codecov](https://codecov.io/gh/mmikhasenko/AlgebraPDF.jl/branch/master/graph/badge.svg)](https://codecov.io/gh/mmikhasenko/AlgebraPDF.jl)

Basic functionality:
 * Attach default values of parameters to a function
 * Update, fix, release parameters
 * constructing a complex model object from set of function:
     - algebra of functions with parameters, e.g. `f₁ + f₂`, or `abs2(f)`.
 * On-fly normalization
 * construction of mixed models in the form `f₁ PDF₁ + f₂ PDF₂ + (1-f₁-f₂) PDF₃`.
 * fitting data distribution using the maximum likelihood (`Optim.jl`)
 * plotting recipes

## Functions with parameters
It is just a function to which a container with parameters (default values) is attached.
The container can be static `NamedTuple`, or extended which can flag parameters as `free` and `fixed`.

There are three main constructors:
1. lambda-function is explicitly given
```julia
FunctionWithParameters(f::F, p::P)
```

2. user-defined stucture which is subtype of `AbstractFunctionWithParameters`:
```julia
struct myAmazingF{P} <: AbstractFunctionWithParameters
    p::P
end
func(d::myAmazingF, x::Number; p=pars(d)) = ... # expression
```

3. using a macro `@makefuntype`:
```julia
@makefuntype myAmazingF(x;p) = ... # expression
```

## Normalized functions

The idea is to attach also the limit to the function and compute the integral of it for the given parameter on-fly.
To make the normalization efficient, a call of the function on the `AbstractVector` implements a broadcasting with a single computation of normalization.
```julia
mypdf(1.1) # use default values of parameters, calls normalization once
mypdf(rand(100)) # use default values of parameters, also calls normalization once
mypdf(rand(100); p = (a=1.2, b=3.3)) # ignors defalt parameters
```

The PDF has two main representations (the ways to define):
1. A struct with the reference to the `unnormdensity<:AbstractFunctionWithParameters`. 
```julia
struct PDFWithParameters{T<:AbstractFunctionWithParameters,L} <: AbstractPDF{1}
    unnormdensity::T
    lims::L
end
```
A regular function can be wrapper to FunctionWithParameters: `FunctionWithParameters((x;p)->p.c0+p.c1*x, (c0=1.0, c1=2.0))`

2. Alternativerly, the density can be defined using a dispatch on a `customaty_type <: AbstractPDF{1}`. E.g.,
```julia
struct Pol1SinSq{T,N} <: AbstractPDF{1}
    p::T
    lims::N
end
func(d::Pol1SinSq, x::Number; p=pars(d)) = p.a*sin(x+p.b)^2+1  # an example of the function
```

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

## Convolution

The most common case of smearing a function with gaussian denisity is implemented.
The convolved function is created with
```julia
f_conv = convGauss(f::F, σ::T) where F <: AbstractFunctionWithParameters
```
σ can be a number, but can also be a function `<: AbstractFunctionWithParameters`.

A customary confolved function or pdf can be defined the same was as e.g. [`FBreitWignerConvGauss`](src/densities.jl).
