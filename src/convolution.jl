
conv_with_gauss(e,f,σ; nσ=5) = conv_f_with_g(e,f,x->standardgauss(x,σ); Δ=nσ*σ)[1]
conv_with_gauss(d::pdf,σ; nσ=5) = pdf((e;p)->conv_with_gauss.(e,x->func(d,x; p=p),σ; nσ=nσ); p=pars(d), lims=lims(d))
# 
conv_f_with_g(x,f,g; Δ::Float64=5.0) = quadgk(y->f(y)*g(x-y), x-Δ, x+Δ)[1]

function conv_with_gauss_sampling(e,f,σ; Ns=10)
    points = range(-3σ, +3σ, length=Ns)
    norm = sum(standardgauss.(points, σ))
    return sum(f(e .+ x) .* standardgauss.(x, σ) for x in points) / norm
end
conv_with_gauss_sampling(d::pdf,σ; Ns=10) =
    pdf((e;p)->conv_with_gauss_sampling(e, x->func(d,x; p=p), σ; Ns=Ns); p=d.p, lims=lims(d))


#

struct convGauss{T,S} <: AbstractPDF
    pdf::T
    σ::S
end

pars(d::convGauss) = pars(d.pdf) + pars(d.σ)
lims(d::convGauss) = lims(d.pdf)
copy(d::convGauss,p) = convGauss(copy(d.pdf,p), d.σ)

function func(d::convGauss, x::Number; p=pars(d))
    σ = func(d.σ, x; p=p)
    g(z) = AlgebraPDF.standardgauss(z, σ)
    f(z) = func(d.pdf, z; p=p)
    return quadgk(y->f(x-y) * g(y), -5*σ, +5*σ)[1]
end