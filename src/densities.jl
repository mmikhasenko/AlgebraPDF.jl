
standardgauss(x,σ) = exp(-x^2/(2*σ^2)) / sqrt(2π) / σ
# 
# shortcut with fixed parameter names
# @makefuntype FGauss(x;p) = exp(-(x-p.μ)^2/(2*p.σ^2))
# 
# better: just order of parameters is fixed
"""
    FGauss{P} <: AbstractFunctionWithParameters

A Gaussian (or normal) function, represented as a subtype of AbstractFunctionWithParameters.
The struct is parametered by the type of its parameters

This struct represents a Gaussian distribution with mean `μ` and standard deviation `σ`. The function is defined as:

```math
f(x) = \\frac{1}{\\sqrt{2π}σ} e^{-\\frac{(x-μ)^2}{2σ^2}}
```

# Fields
- `p::P`: Parameters of the Gaussian distribution. `P` is a type, a named tuple in the simplest case,
that should contain the keys for the mean and standard deviation.

# Example
```julia
gaussian = FGauss((my_μ=0.0, my_σ=1.0))
gaussian(1.1)  # Evaluate the Gaussian function at x=1.1
gaussian(1.1; p=(my_μ=0.5, my_σ=1.5))  # Evaluate with adjusted parameters.
```

# Notes
- The order of parameters in `P` is important: first mean, then standard deviation.
- To evaluate the function at a given point `x`, one can do either `gaussian(x)`, or `func(gaussian, x)`.
"""
struct FGauss{P} <: AbstractFunctionWithParameters
    p::P
end
function func(d::FGauss, x::NumberOrTuple; p=pars(d))
    μ,σ = (getproperty(p,s) for s in keys(d.p))
    standardgauss(x-μ,σ)
end

###############################################################

amplitudeBW(x,m,Γ) = m*Γ/(m^2-x^2-1im*m*Γ)
amplitudeBWsq(x,m,Γ) = abs2(amplitudeBW(x,m,Γ))

# shortcut with fixed parameter names
# @makefuntype FBreitWigner(x;p) = amplitudeBW(x,p.m,p.Γ)
# 
# better: just order of parameters is fixed
struct FBreitWigner{P} <: AbstractFunctionWithParameters
    p::P
end
function func(d::FBreitWigner, x::NumberOrTuple; p=pars(d))
    m,Γ = (getproperty(p,s) for s in keys(d.p))
    amplitudeBW(x,m,Γ)
end

###############################################################

# shortcut with fixed parameter names
# @makefuntype FExp(x;p) = exp(-(x-p.μ)^2/(2*p.σ^2))
# 
# better: just order of parameters is fixed

"""
    FExp{P} <: AbstractFunctionWithParameters

An Exponential function, represented as a subtype of AbstractFunctionWithParameters.
The struct is parameterized by the type of its parameters.

This struct represents an exponential distribution function with parameter `β`. The function is defined as:

```math
f(x) = e^{βx}
```

# Fields
- `p::P`: Parameters of the exponential function. `P` is a type, typically a named tuple, that should contain an entry for the slope.

# Example
```julia
expfunc = FExp((my_β=-0.5))
expfunc(1.1)  # Evaluate the Exponential function at x=1.1
expfunc(1.1; p=(my_β=-0.7))  # Evaluate with an adjusted parameter.
```

# Notes
- The parameter `β` defines the rate or scale of the exponential function.
- To evaluate the function at a given point `x`, one can use either `expfunc(x)` or `func(expfunc, x)`.
"""
struct FExp{P} <: AbstractFunctionWithParameters
    p::P
end
function func(d::FExp, x::NumberOrTuple; p=pars(d))
    β, = (getproperty(p,s) for s in keys(d.p))
    exp(x*β)
end

###############################################################

# shortcut with fixed parameter names
# @makefuntype FPowExp(x;p) = x^p.α*exp(x*p.β)
# 
# better: just order of parameters is fixed
struct FPowExp{P} <: AbstractFunctionWithParameters
    p::P
end
function func(d::FPowExp, x::NumberOrTuple; p=pars(d))
    α,β = (getproperty(p,s) for s in keys(d.p))
    x^α*exp(x*β)
end

###############################################################

# shortcut with fixed parameter names
# @makefuntype FPol1(x;p) = p.c0+x*p.c1
# @makefuntype FPol2(x;p) = p.c0+x*p.c1+x^2*p.c2
# 
# better: just order of parameters is fixed
"""
    FPol{P} <: AbstractFunctionWithParameters

A polynomial function, a linear sum of power serience in `x` with an increamental power starting from 0. 
The coefficients are defined in the parameters. The function is generally defined as:

```math
f(x) = c_0 + c_1 x + c_2 x^2 + c_3 x^3 + \\cdots
```

# Fields
- `p::P`: Parameters of the polynomial function.
`P` is a type, typically a named tuple, that should contain the coefficients of the polynomial, e.g. `c0`, `c1`, `c2`.

# Example
```julia
polynomial = FPol((c0=1.0, c1=-0.5, c2=0.25))
polynomial(2.0)  # Evaluate the polynomial function at x=2.0
polynomial(2.0; p=(c0=1.0, c1=-0.5, c2=0.1))  # Evaluate with adjusted function
```

# Notes
- The coefficients in `P` follow the order of the polynomial terms, starting from the constant term (c0) to higher-degree terms (c1, c2, ...).
The user is reponsible to use appropriate names for the coeffients. `FPol((c6=1.0, c2=-0.5, c0=0.25))` is a valid construction, however, extremely comfusing as referres to `c_8 + c_2 x + c_0 x^2` 
- To evaluate the polynomial function at a given point `x`, one can use either `polynomial(x)` or `func(polynomial, x)`.
"""
struct FPol{P} <: AbstractFunctionWithParameters
    p::P
end
function func(d::FPol, x::NumberOrTuple; p=pars(d))
    cs = (getproperty(p,s) for s in keys(d.p))
    sum(x^(i-1)*c for (i,c) in enumerate(cs))
end

###############################################################
standarddoublegauss(x,σ,r,n) =
    r*standardgauss(x,σ) + (1-r)*standardgauss(x,n*σ)

# shortcut with fixed parameter names
# @makefuntype FDoubleGaussFixedRatio(x;p) = standarddoublegauss(x-p.μ,p.σ,p.r,p.n)
# 
# better: just order of parameters is fixed
"""
    FDoubleGaussFixedRatio{P} <: AbstractFunctionWithParameters

A combination of two Gaussian functions with a fixed ratio, represented as a subtype of AbstractFunctionWithParameters. 
The struct is parameterized by the type of its parameters, which define the characteristics of the two Gaussian functions.

This struct represents a function that is a weighted sum of two Gaussian distributions with a common mean `μ` but different standard deviations `σ` and `n*σ`. The ratio `r` defines the weight of the first Gaussian in the sum. The function is defined as:

```math
f(x) = r \\frac{1}{\\sqrt{2π}σ} e^{-\\frac{(x-μ)^2}{2σ^2}} + (1-r) \\frac{1}{\\sqrt{2π}(nσ)} e^{-\\frac{(x-μ)^2}{2(nσ)^2}}
```

# Fields
- `p::P`: Parameters of the double Gaussian function.
`P` is a type, typically a named tuple, that should contain the keys for mean, standard deviation of the first Gaussian,
ratio of the first Gaussian, and the ratio of the second Gaussian's standard deviation to the first' one.

# Example
```julia
doubleGauss = FDoubleGaussFixedRatio((μ=0.0, σ=1.0, r=0.5, n=2.0))
doubleGauss(1.1)  # Evaluate the double Gaussian function at x=1.1
doubleGauss(1.1; p=(μ=0.5, σ=1.5, r=0.7, n=2.5))  # Evaluate with adjusted parameters.
```
"""
struct FDoubleGaussFixedRatio{P} <: AbstractFunctionWithParameters
    p::P
end
function func(d::FDoubleGaussFixedRatio, x::NumberOrTuple; p=pars(d))
    μ,σ,r,n = (getproperty(p,s) for s in keys(d.p))
    standarddoublegauss(x-μ,σ,r,n)
end

###############################################################

struct FBreitWignerConvGauss{P} <: AbstractFunctionWithParameters
    p::P
end
function func(d::FBreitWignerConvGauss, x::NumberOrTuple; p=pars(d))
    m,Γ,σ = (getproperty(p,s) for s in keys(d.p))
    # 
    quadgk(
        y->abs2(amplitudeBW(x-y, m, Γ)) *
            standardgauss(y,σ), -5*σ, +5*σ)[1]
end

###############################################################

struct FTabulated{T} <: AbstractFunctionWithParameters
    itr::T
end
FTabulated(xv,yv) = FTabulated(interpolate((xv,), yv, Gridded(Linear())))
# 
function func(d::FTabulated, x::NumberOrTuple; p=pars(d))    
    tlims = extrema(collect(d.itr.knots)[1])
    inrange(x,tlims) ? d.itr(x) : 0.0
end
pars(d::FTabulated, isfree::Bool) = ∅

###############################################################
###############################################################

abstract type AbstractCrystalBall <: AbstractFunctionWithParameters end
struct FLeftSideCrystalBall{P} <: AbstractCrystalBall
    p::P
end
function func(d::FLeftSideCrystalBall, x; p=pars(d))
    μS, σS, αS, nS = keys(d.p)
    μ, σ, α, n = getproperty.(Ref(d.p), (μS, σS, αS, nS))
    # 
    absα = abs(α)
    C = n/absα / (n-1) * exp(-absα^2/2)
    D = sqrt(π/2)*(1+erf(absα/sqrt(2)))
    N = 1/(σ*(C+D))
    # 
    x̂ = (x-μ) / σ
    x̂ > -α && return exp(-x̂^2/2) * N
    # 
    A = (n/absα)^(n) * exp(-absα^2/2)
    B = n/absα-absα
    return A*(B-x̂)^(-n) * N
end

struct FRightSideCrystalBall{P} <: AbstractCrystalBall
    p::P
end
function func(d::FRightSideCrystalBall, x; p=pars(d))
    μS, σS, αS, nS = keys(d.p)
    μ, σ, α, n = getproperty.(Ref(d.p), (μS, σS, αS, nS))
    # 
    absα = abs(α)
    C = n/absα / (n-1) * exp(-absα^2/2)
    D = sqrt(π/2)*(1+erf(absα/sqrt(2)))
    N = 1/(σ*(C+D))
    # 
    x̂ = (x-μ) / σ
    x̂ < α && return exp(-x̂^2/2) * N
    # 
    A = (n/absα)^(n) * exp(-absα^2/2)
    B = n/absα-absα
    return A*(B+x̂)^(-n) * N
end


struct FDoubleSideCrystalBall{P} <: AbstractCrystalBall
    p::P
end
function func(d::FDoubleSideCrystalBall, x; p=pars(d))
    μS, σS, αS, nS = keys(d.p)
    μ, σ, α, n = getproperty.(Ref(d.p), (μS, σS, αS, nS))
    # 
    absα = abs(α)
    C = n/absα / (n-1) * exp(-absα^2/2)
    D = sqrt(π/2)*erf(absα/sqrt(2))
    N = 1/(σ*2*(C+D))
    # 
    x̂ = (x-μ) / σ
    abs(x̂) < α && return exp(-x̂^2/2) * N
    # 
    A = (n/absα)^(n) * exp(-absα^2/2)
    B = n/absα-absα
    return A*(B+abs(x̂))^(-n) * N
end

import Base:show
function show(io::IO, d::AbstractCrystalBall)
    μS, σS, αS, nS = keys(d.p)
    μ, σ, α, n = getproperty.(Ref(d.p), (μS, σS, αS, nS))
    # 
    println(io, "$(typeof(d)[1])((")
    println(io, "    $(μS)=$(μ), # mode")
    println(io, "    $(σS)=$(σ), # sigma")
    println(io, "    $(αS)=$(α), # alpha")
    println(io, "    $(nS)=$(n))) # n")
end