
standardgauss(x,σ) = exp(-x^2/(2*σ^2)) / sqrt(2π) / σ
amplitudeBW(x,m,Γ) = m*Γ/(m^2-x^2-1im*m*Γ)
amplitudeBWsq(x,m,Γ) = abs2(amplitudeBW(x,m,Γ))
# 
function aGauss(p, lims)
    μ, σ = keys(p)
    return pdf((x;p)->standardgauss.(x .- getproperty(p,μ), getproperty(p,σ)), p, lims)
end
function aBreitWigner(p, lims)
    m, Γ = keys(p)
    return pdf((x;p)->amplitudeBWsq.(x,getproperty(p,m), getproperty(p,Γ)), p, lims)
end
function aExp(p, lims)
    α, = keys(p)
    return pdf((x;p)->exp.(x .*getproperty(p,α)), p, lims)
end
function aPowExp(p, lims)
    α,β = keys(p)
    return pdf((x;p)->x.^getproperty(p,α) .* exp.(x .*getproperty(p,β)), p, lims)
end
function aPol(p, lims)
    cs = keys(p)
    return pdf((x;p)->sum(x.^(i-1) .* getproperty(p,c) for (i,c) in enumerate(cs)), p, lims)
end

standarddoublegauss(x,σ,r,n) =
    r*AlgebraPDF.standardgauss(x,σ) + (1-r)*AlgebraPDF.standardgauss(x,n*σ)

function aDoubleGaussFixedRatio(pars, lims; fixpars)
    μ,σ = keys(pars)
    r,n = fixpars
    return pdf((x;p)->standarddoublegauss.(x .- getproperty(p, μ), getproperty(p, σ), r,n),
        pars, lims)
end

function aBreitWignerConvGauss(pars, lims; fixpars)
    m,Γ = keys(pars)
    σ, = fixpars
    density = (x;p)->conv_with_gauss.(x,
        y->abs2(AlgebraPDF.amplitudeBWsq(y, getproperty(p,m), getproperty(p,Γ))), σ)
    return pdf(density, pars, lims)
end

function aTabulated(xv,yv,lims)
    itr = interpolate((xv,), yv, Gridded(Linear()))
    tlims = (prevfloat(xv[1]),nextfloat(xv[end]))
    f(x) = inrange(x,tlims) ? itr(x) : 0.0
    fixedshapepdf(x->f.(x), lims)
end
