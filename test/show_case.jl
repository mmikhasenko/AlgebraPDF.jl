using Test
using AlgebraPDF

snl = pdf(@. (x;p) -> exp(-(x-p.μ)^2/(2*p.σ^2)); p0 = (μ=1.4, σ=0.15), lims=(0, 3))
bkg = pdf(@. (x;p) -> sqrt(x)*exp(-p.α*x); p0 = (α=1.3,), lims=(0, 3))
bkg *= (fb=2.5,)
pdf_sum = snl + bkg

# generating
const data = generate(10000, pdf_sum; p=pdf_sum.p0);
@test length(data) == 10000

# fitting
pfr = v2p(fit_llh(data, pdf_sum;
  init_pars=p2v(pdf_sum.p0, pdf_sum)), pdf_sum)
@test length(pfr) == 4
