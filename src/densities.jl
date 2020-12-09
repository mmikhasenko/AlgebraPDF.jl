
standardgauss(x,σ) = exp(-x^2/(2*σ^2)) / sqrt(2π) / σ
standardBW(x,m,Γ) = abs2(m*Γ/(m^2-x^2-1im*m*Γ))
# 
function aGauss(p, lims)
    μ, σ = keys(p)
    return pdf((x;p)->standardgauss.(x .- getproperty(p,μ), getproperty(p,σ)), p, lims)
end
function aBreitWigner(p, lims)
    m, Γ = keys(p)
    return pdf((x;p)->standardBW.(x,getproperty(p,m), getproperty(p,Γ)), p, lims)
end
function aExp(p, lims)
    α, = keys(p)
    return pdf((x;p)->exp.(x .*getproperty(p,α)), p, lims)
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
        y->abs2(AlgebraPDF.standardBW(y, getproperty(p,m), getproperty(p,Γ))), σ)
    return pdf(density, pars, lims)
end

function aTabulated(xv,yv,lims)
    itr = interpolate((xv,), yv, Gridded(Linear()))
    fixedshapepdf(x->itr.(x), lims)
end
