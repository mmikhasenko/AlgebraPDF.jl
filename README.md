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


## How to call
The object behave similar to a regular function with a keyword argument `p` set to `freepars(d)` by default.
Once `p` is used a full set of parameters needs to be provided. 
```julia
g = FGauss((μ=1.2, σ=0.1))
# call on a single argument
g(0.7) # use default parameters
g(0.7; p=(μ=1.6, σ=0.2)) # use provided parameters
g(0.7, [1.6, 0.2]) # same as before, see `NamedTuple{keys(freepars(g))}`
# incorrect calls
g(0.7; p=(μ=1.6,)) # error: no σ is given
g(0.7; p=(σ=0.2,)) # error: no μ is given
# broadcasting
g(rand(10))
g(rand(10); p=(μ=1.6, σ=0.2))
```

## Update parameter values
The function structure immutable, therefore, the method on modification of parameters returns a new object
```julia
d = FunctionWithParameters((x;p)->x^2+p.a, (a=1.0,))
# 
d′ = updatepar(d, :a, 2.0)
d′′ = updatepars(d, (a=2.0,))
```

## Fix/release parameter values

By extending the parameter structure, one gets a possibility to fix/release the parameters.
```julia
g = FGauss(Ext(μ=1.2, σ=0.1)) # Ext for extended
# 
g′ = fixpar(g, :μ)
g′ = fixpar(g, :μ, 1.3)
# 
freepars(g′) # (σ=0.1,)
fixedpars(g′) # (μ=1.3,)
#
g′(0.7; p=(σ=0.2,)) # now works since μ is fixed
g′(0.7; p=(σ=0.2, μ=1.6)) # will use the fixed value of μ=1.2
g′′ = updatepar(g′, :μ, 1.2) # however, the update fill do
# 
g′′′ = releasepar(g′′, :μ)
g′′′ == g # true
```

## Fit free parameters
The function can be subjected to the unbinned fit.
```julia
fit_results = fit_llh(data, f; init)
@unpack best = fit_results
```
The method returns a named tuple with the named tuple of optimized parameters,
parameters with uncertainties, and best-estimation function.

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

## Fixed parameter names or fixed parameter order

Creating a function or pdf can be conveniently done with macro
```julia
# for BW1 <: AbstractPDF{1}
@makepdftype BW1(x, p) = p.m*p.Γ/(p.m^2-x^2-1im*p.m*p.Γ)

# for BW1 <: AbstractPDF{1}
@makefuntype BW1(x, p) = p.m*p.Γ/(p.m^2-x^2-1im*p.m*p.Γ)
```
The latter expands into
```julia
# implementation with NAMES of parameters build into the funciton call
struct BW1{P} <: AbstractFunctionWithParameters
    p::P
end
import AlgebraPDF:func
func(bw::BW1, x::Number; p=pars(bw)) = p.m*p.Γ/(p.m^2-x^2-1im*p.m*p.Γ)
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
