
standardgauss(x,σ) = exp(-x^2/(2*σ^2)) / sqrt(2π) / σ
# 
# shortcut with fixed parameter names
# @makefuntype FGauss(x;p) = exp(-(x-p.μ)^2/(2*p.σ^2))
# 
# better: just order of parameters is fixed
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