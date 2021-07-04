
standardgauss(x,σ) = exp(-x^2/(2*σ^2)) / sqrt(2π) / σ
# 
# shortcut with fixed parameter names
# @makefuntype FGauss(x;p) = exp(-(x-p.μ)^2/(2*p.σ^2))
# 
# better: just order of parameters is fixed
struct FGauss{P} <: AbstractFunctionWithParameters
    p::P
end
function func(d::FGauss, x::Number; p=pars(d))
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
function func(d::FBreitWigner, x::Number; p=pars(d))
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
function func(d::FExp, x::Number; p=pars(d))
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
function func(d::FPowExp, x::Number; p=pars(d))
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
function func(d::FPol, x::Number; p=pars(d))
    cs = (getproperty(p,s) for s in keys(d.p))
    sum(x^(i-1)*c for (i,c) in enumerate(cs))
end

###############################################################

standarddoublegauss(x,σ,r,n) =
    r*standardgauss(x,σ) + (1-r)*standardgauss(x,n*σ)

# shortcut with fixed parameter names
# @makefuntype FDoubleGaussFixedRatio(x;p) = gauss(x,p.σ,p.r,p.n)
# 
# better: just order of parameters is fixed
struct FDoubleGaussFixedRatio{P} <: AbstractFunctionWithParameters
    p::P
end
function func(d::FDoubleGaussFixedRatio, x::Number; p=pars(d))
    μ,σ,r,n = (getproperty(p,s) for s in keys(d.p))
    standarddoublegauss(x-μ,σ,r,n)
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
function func(d::FDoubleGaussFixedRatio, x::Number; p=pars(d))
    μ,σ,r,n = (getproperty(p,s) for s in keys(d.p))
    standarddoublegauss(x-μ,σ,r,n)
end

###############################################################

struct FBreitWignerConvGauss{P} <: AbstractFunctionWithParameters
    p::P
end
function func(d::FBreitWignerConvGauss, x::Number; p=pars(d))
    m,Γ,σ = (getproperty(p,s) for s in keys(d.p))
    conv_with_gauss.(x, y->abs2(amplitudeBW(y, m, Γ)), σ)
end

###############################################################

struct FTabulated{T} <: AbstractFunctionWithParameters
    itr::T
end
FTabulated(xv,yv) = FTabulated(interpolate((xv,), yv, Gridded(Linear())))
# 
function func(d::FTabulated, x::Number; p=pars(d))    
    tlims = extrema(collect(d.itr.knots)[1])
    inrange(x,tlims) ? d.itr(x) : 0.0
end
pars(d::FTabulated, isfree::Bool) = ∅

###############################################################
