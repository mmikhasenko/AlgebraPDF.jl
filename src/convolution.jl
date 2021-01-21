
conv_with_gauss(e,f,σ; nσ=3) = quadgk(y->f(y)*standardgauss(e-y,σ), e-nσ*σ, e+nσ*σ)[1]
conv_with_gauss(d::pdf,σ; nσ=3) = pdf((e;p)->conv_with_gauss.(e,x->func(d,x; p=p),σ; nσ=nσ); p=pars(d), lims=lims(d))

function conv_with_gauss_sampling(e,f,σ; Ns=10)
    points = range(-3σ, +3σ, length=Ns)
    norm = sum(standardgauss.(points, σ))
    return sum(f(e .+ x) .* standardgauss.(x, σ) for x in points) / norm
end
conv_with_gauss_sampling(d::pdf,σ; Ns=10) =
    pdf((e;p)->conv_with_gauss_sampling(e, x->func(d,x; p=p), σ; Ns=Ns); p=d.p, lims=lims(d))
