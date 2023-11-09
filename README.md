# AlgebraPDF


| **Documentation**                                                               | **Build Status**                                                                                |
|:-------------------------------------------------------------------------------:|:-----------------------------------------------------------------------------------------------:|
| [![][docs-stable-img]][docs-stable-url] [![][docs-dev-img]][docs-dev-url] | [![][GHA-img]][GHA-url] [![][codecov-img]][codecov-url] |



[docs-dev-img]: https://img.shields.io/badge/docs-dev-blue.svg
[docs-stable-img]: https://img.shields.io/badge/docs-stable-blue.svg
[GHA-img]: https://github.com/JuliaDocs/Documenter.jl/workflows/CI/badge.svg
[codecov-img]: https://codecov.io/gh/mmikhasenko/AlgebraPDF.jl/branch/master/graph/badge.svg


[docs-dev-url]: https://mmikhasenko.github.io/AlgebraPDF.jl/dev/
[docs-stable-url]: https://mmikhasenko.github.io/AlgebraPDF.jl/stable/
[GHA-url]: https://github.com/mmikhasenko/AlgebraPDF.jl/actions?query=workflows/CI
[codecov-url]: https://codecov.io/gh/mmikhasenko/AlgebraPDF.jl



## Basic functionality:
 * Attach default values of parameters to a function
 * Update, fix, release parameters
 * constructing a complex model object from set of function:
     - algebra of functions with parameters, e.g. `f₁ + f₂`, `abs2(f)`, or `log(f)`.
 * On-fly normalization
 * construction of mixed models in the form `f₁ PDF₁ + f₂ PDF₂ + f₃ PDF₃`.
 * construction of likelihood function and extended likelihood function
 * plotting recipes

Current implementation is limited to immutable operations.

## An example 
![plots/sixcompfit.png](plots/sixcompfit.png)
### Model construction and Extended MLE fit
```julia
using AlgebraPDF
using DelimitedFiles


const xth = 2.95
const support = xth .+ (0, 0.22)

data = readdlm("data.txt")[:,1]
# Nev = size(data,1)

# amplitude for the signal
Φ2(x) = sqrt(x-xth)
Γ(x,m,Γ₀) = Γ₀*Φ2(x)/Φ2(m)
breitwigner(x,m,Γ₀) = 1/(m^2-x^2-1im*m*Γ(x,m,Γ₀)) 

# I) phase space function, also the background
phasespace = FunctionWithParameters((x;p)->Φ2(x), ∅) # pass λ-function
backgrpdf = Normalized(phasespace, support) # get PDF

# II) define a type SimpleBW and the method `func` for dispatch
struct SimpleBW{P} <: AbstractFunctionWithParameters
    p::P
end
import AlgebraPDF:func
function func(bw::SimpleBW, x::NumberOrTuple; p=pars(bw))
    m,Γ = (getproperty(p,s) for s in keys(bw.p))
    breitwigner(x, m, Γ)
end

# Signal1-4: |A|^2 * phase_space
signalpdfs = 
    [Normalized(abs2(A)*phasespace, support)
        for A in [
            SimpleBW((m1=3.00, Γ1=6.5e-3)),
            SimpleBW((m2=3.05, Γ2=2.3e-3)),
            SimpleBW((m3=3.06, Γ3=4.0e-3)),
            SimpleBW((m4=3.09, Γ4=9.9e-3))]
        ]

# III) The threshold function, also |Ath|^2 * phase_space
@makefuntype SimpleBWg(x;p) =
    1/(p.m0^2 - x^2 - 1im*p.g^2*Φ2(x))
#
signalpdf0 = Normalized(
        abs2(SimpleBWg((m0=2.95, g=0.01)))*phasespace,
        support)

# the full model - sum of components
model0 = 
    signalpdf0 * (f0=0.05Nev,) +
    signalpdfs[1] * (f1=0.25Nev,) + signalpdfs[2] * (f2=0.15Nev,) + 
    signalpdfs[3] * (f3=0.15Nev,) + signalpdfs[4] * (f4=0.55Nev,) +
    backgrpdf * Ext(fb=0.3Nev,) # Ext - to be able to fix
# 
model1 = fixpar(model0, :fb, 17.2)
# 
# fit with your favorite package
```
The plotting commands see in [plots/example.jl](plots/example.jl). Detailed description of the methods follows.

## Call the function
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

## Algebra of functions

Several operations on a function are implemented, as `abs2` and `log`.
```julia
# an example of complex-valued function
f = FunctionWithParameters((x;p)->1/(p.m^2-x^2-1im*p.m*p.Γ), (m=0.77,Γ=0.15))

fabs2 = abs2(f1)
flog = log(fabs2)
```

A simple arithmetics on a pair of function also works.
```julia
# operations with two functions
f1 = FunctionWithParameters((x;p)->abs(x)<1 ? 1-(2x)^2 : -3, ∅)
f2 = FunctionWithParameters((x;p)->p.a*exp(-x)+p.b, (a=0.1,b=0.1))
#
fprod = f1*f2
# 
fsum = f1+f2
fsub = f1-f2
```
The summation and subtraction adds a new parameter for the coefficient of the functions, `p.α1*f1+p.α2*f2`.
The parameter can be passed with 
```julia
+(f1,f2; p=(c1=1.1,c2=2.2))
+(f1,f2; p=Ext(c1=1.1,c2=2.2)) # the coefficients can be fixed 
# using explicit constructor:
FSum([f1,f2], (c1=1.1,c2=2.2))
FSum([f1,f2], Ext(c1=1.1,c2=2.2))
# subtraction is an addition
f1-f2 == +(f1,f2; p=(α1=1.0,α2=-1.0)) # true
``` 
Perhaps, a more transparent constuction of the same sum can be done using a linear decomposition:
```julia
(c1=1.1,) * f1 + (c2=2.2,) * f2
```
Where a multiplication of a function to a parameter type returns a sum with a single term.
A special method on the addition on the sum is called.

## Create/implement the functions with parameters
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
func(d::myAmazingF, x::NumberOrTuple; p=pars(d)) = ... # expression
```

3. using a macro `@makefuntype`:
```julia
@makefuntype myAmazingF(x;p) = ... # expression
```

## Create/implement normalized functions

The idea is to attach also the limit to the function and compute the integral of it for the given parameter on-fly.
To make the normalization efficient, a call of the function on the `AbstractVector` implements a broadcasting with a single computation of normalization.
```julia
myNormalized(1.1) # use default values of parameters, calls normalization once
myNormalized(rand(100)) # use default values of parameters, also calls normalization once
myNormalized(rand(100); p = (a=1.2, b=3.3)) # ignors defalt parameters
```

The PDF has two main representations (the ways to define):
1. A struct with the reference to the `unnormdensity<:AbstractFunctionWithParameters`. 
```julia
struct Normalized{T<:AbstractFunctionWithParameters,L} <: AbstractPDF{1}
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
func(d::Pol1SinSq, x::NumberOrTuple; p=pars(d)) = p.a*sin(x+p.b)^2+1  # an example of the function
```
The limits can be checked with `lims(d)`.

## Defined parameter names or defined parameter order

Creating a function or pdf can be conveniently done with macro
```julia
# for BW1 <: AbstractPDF{1}
@makepdftype BW1(x, p) = p.m*p.Γ/(p.m^2-x^2-1im*p.m*p.Γ)
```
which expands into
```julia
# implementation with NAMES of parameters build into the funciton call
struct BW1{P} <: AbstractFunctionWithParameters
    p::P
end
import AlgebraPDF:func
func(bw::BW1, x::NumberOrTuple; p=pars(bw)) = p.m*p.Γ/(p.m^2-x^2-1im*p.m*p.Γ)
```

Slightly better implementation where only the order of the arguments are fixed, while the names are determined when the instance is created.
```julia
# implementation with ORDER of parameters build into the funciton call
struct BW2{P} <: AbstractFunctionWithParameters
    p::P
end
import AlgebraPDF:func
function func(bw::BW2, x::NumberOrTuple; p=pars(bw))
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

## Some defined functions

For convenience, some standard functions are predefined.
```julia
FGauss((μ=1.2, σ=0.2))
FDoubleGaussFixedRatio((μ=1.2, σ=0.2, r=0.2, n=5))
FBreitWigner((m=0.77, Γ=0.15))
FExp((τ=-0.1,))
FPowExp((n=1.2,τ=-0.1))
FPol((a0=1, a1=2, a2=3, a3=5))
FBreitWignerConvGauss((μ=0.77, Γ=0.15, σ=2.2))
# 
xv = -π/2:0.1:π/2
yv = sin.(xv)
FTabulated(xv,yv)
```

The corresponding pdf can be defined with `Normalized` by adding a limits. E.g.,
```julia
d = Normalized(FGauss((μ=1.2, σ=0.2)), (-1,4))
```

## Higher dimensions

Function with higher dimensions expect the variable provided as a `Tuple`, i.e. `(x,y)`, otherwise,
the construction and usage are analogous to the one dimensional case.
```julia
@makefuntype Amazing2D(x;p) = (x[1]-p.x0)^2+(x[2]-p.y0)^2-p.R0^2
a = Amazing2D((x0=1.1, y0=2.1, R0=0.0))
a((1.1,2.1)) # returns 0
# 
data = collect(zip(rand(10), rand(10)))
a(data) # return a vector of 10 elements

@makefuntype Amazing3D(x;p) = (x[1]-p.x0)^2+(x[2]-p.y0)^2+(x[3]-p.z0)^2-p.R0^2
b = Amazing3D((x0=1.1, y0=2.1, z0=3.1, R0=-1.0))
b((1.1,2.1,3.1)) == -1
# 
data3d = collect(zip(rand(10), rand(10), rand(10)))
b(data3d; p=(x0=1, y0=2, z0=3, R0=-1) )  # return a vector of 10 elements
```
## Plotting
The function object can be plotted as a regular function,
```julia
f1 = FGauss((μ1=1.2, σ1=0.2))
plot(f1, 0, 5)
```
It is replaced to a lambda function by the type recipe.

For a PDF that has the limits, the plotting command will just work
```julia
d1 = Normalized(FGauss((μ1=1.2, σ1=0.2)), (0,4))
plot(d1, l=(:orange,3),
    lab="FGauss(μ1=1.2, σ1=0.2)",
    title="Gaussian normalized in (0,4)")
```
The `normalization` value other than unit can be passed with the plotting method
```julia
plot(d1, normalization=1.0, Nsample=100)
```
where `Nsample` is the number of points at which the function is sampled.

Plotting recipe for 2d functions is defined.
```julia
@makefuntype Amazing2D(x;p) = (x[1]-p.x0)^2+(x[2]-p.y0)^2-p.R0^2
a = Amazing2D((x0=0.1, y0=-0.1, R0=0.0))
xv, yv = -1:0.1:2, -2:0.1:1
heatmap(xv, yv, updatepar(a, :x0, -0.5))
contour(xv, yv, log(abs2(a)))
surface(xv, yv, abs2(a))
```

## Sample data from PDF
A data set distributed according a given model (`<:AbstractFunctionWithParameter`, and defined limits)
can be generated using
a function `generate(model, N)` with `N` being the required number of events
```julia
@makepdftype ExpAbs(x;p) = exp(-abs(x)/p.σ)
model = ExpAbs((σ=2.1,), (-3,2))
# 
data1 = generate(model, 1000)
data2 = generate(model, 1000; p=(σ=0.3,), Nbins=1000)
```
The inversion of the cdf is done numerically using the binning controlled by `Nbins` variable.
A set of parameters to be used can be passed as a key argument.
The function `rand` will also work:
```julia
value = rand(model) # not efficient
data3 = rand(model, 1000) # same as generate(model, 1000)
```
however, a call for a single number in inefficient (the integral is not stored but computed every time),
the call with a set size does the same `generate` and does not take key arguments.
Currently, `generate` works only in one dimension.
