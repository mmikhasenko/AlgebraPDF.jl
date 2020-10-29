using AlgebraPDF
using Plots
using TypedTables
using StaticArrays
theme(:wong2)
#                  
g1 = pdf(@. (x;p)->1/p.σ1*exp(-(x-p.μ1)^2/(2*p.σ1^2)); p0=(μ1= 2.1, σ1=0.7 ), lims=(-3, 3))
g2 = pdf(@. (x;p)->1/p.σ2*exp(-(x-p.μ2)^2/(2*p.σ2^2)); p0=(μ2=-0.7, σ2=0.7 ), lims=(-3, 3))
g3 = pdf(@. (x;p)->1/p.σ3*exp(-(x-p.μ3)^2/(2*p.σ3^2)); p0=(μ3= 0.3, σ3=0.1), lims=(-3, 3))
# 
mm0 = MixedModel(SVector(g1, g2, g3), (f1=0.33,f2=0.6))

plot(x->mm0(x),-3,3)
plot!(x->0.33*mm0.components[1](x), -3, 3)
plot!(x->0.60*mm0.components[2](x), -3, 3)
plot!(x->0.07*mm0.components[3](x), -3, 3)

sample = vcat((0.5.*randn(1000) .- 1.0), (0.7.*randn(100) .+ 2.0), 0.2 .* randn(300) .+ 0.3)
sample = filter(sample) do x
    -3<x<3
end
histogram(sample, bins=80)
p0 = (μ1=-1.0, σ1=0.5, μ2=2.1, σ2=0.7, μ3=0.3, σ3=0.2, f1=10/14, f2=1/14)

collectpars(mm0)
fr = fit_llh(sample, mm0, init_pars=p2v(p0, mm0))
pfr = v2p(fr.minimizer, mm0)

stephist(sample,bins=70, norm=true)
plot!(x->mm0(x; p=p0),-3,3)
plot!(x->mm0(x; p=pfr),-3,3)
