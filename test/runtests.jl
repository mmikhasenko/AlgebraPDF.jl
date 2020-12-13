using AlgebraPDF
using StaticArrays
using Test
using SpecialFunctions

@testset "Operations with parameters" begin
    p = (μ = 1.1, σ = 2.2, f = 2.2)

    sp = subtractpars(p, (:μ, :σ))
    @test keys(sp) == (:f,)
    @test sp == subtractpars(p, [:μ, :σ])
    #
    sp = selectpars(p, (:μ, :σ))
    @test keys(sp) == (:μ, :σ)
    @test sp == selectpars(p, [:μ, :σ])
    #
    up = updatepars(  p, (μ = 3.1, σ = 5.2))
    @test length(up) == 3
    @test up.μ == 3.1
    @test up.σ == 5.2
    @test up.f == 2.2
end

@testset "Basic operations" begin
    BW(s, m, Γ) = 1 / (m^2 - s - 1im*m*Γ)
    #
    pdf1 = pdf(@. (e;p)->abs2(BW(e^2, p.m1, p.Γ1));
        p = (m1=0.25, Γ1=2e-3), lims = (0, 0.15))
    #
    @test pdf1(rand()) != 0.0
    @test length(pdf1(rand(10))) == 10
    #
    pdf2 = pdf(@. (e;p)->abs2(BW(e^2, p.m2, p.Γ2));
        p = (m2=0.1, Γ2=14e-3), lims = (0, 0.15))
    #
    x0 = 1.1; v0 = rand(2)
    @test pdf2(x0, v0) == pdf2(x0; p=v2p(v0, pdf2))
    #
    pdf2 *= (f2=3.0,)
    @test length(collectpars(pdf2)) == 3
    #
    pdf_sum = pdf1 + pdf2
    @test npars(pdf_sum) == 5

    pdf_ratio = pdf1 / pdf2
    @test npars(pdf_ratio) == 5
end

function df(f1,f2,lims; Ns = 10)
    xv = range(lims..., length=Ns)
    dyv = f1.(xv) .- f2.(xv)
    return sqrt(sum(abs2, dyv) / Ns)
end

@testset "convolution with gauss" begin

    @test AlgebraPDF.standardgauss(0,1) ≈ 1/sqrt(2π)
    @test AlgebraPDF.standardgauss(2,2) ≈ exp(-1/2)/sqrt(2π*4)
    # tests
    σ0 = 0.3
    aconv(e) = (erf(e/sqrt(2*σ0^2))+1)/2
    #
    nconv(e) = conv_with_gauss(e, x->x>0, σ0)
    @test df(nconv, aconv, (-1, 1); Ns=1000) < 0.01
    #
    step = pdf(@. (e;p)->e>0; p=∅, lims=(-1,1))
    smeared_step = conv_with_gauss(step, σ0)
    nconv(e) = func(smeared_step, e)
    @test df(nconv, aconv, (-1, 1); Ns=1000) < 0.01
    # 
    smeared_step_sampling = conv_with_gauss_sampling(step, σ0; Ns=50)
    nconv(e) = func(smeared_step_sampling, e)
    @test df(nconv, aconv, (-1, 1); Ns=1000) < 0.01
    # 
    smeared_step_sampling = conv_with_gauss_sampling(step, σ0; Ns=10)
    nconv(e) = func(smeared_step_sampling, e)
    @test 0.01 < df(nconv, aconv, (-1, 1); Ns=1000) < 0.04
end

@testset "fix parameters example" begin
    d0 = pdf(@. (e;p)->e^2+p.a; p=(a=1.0,), lims=(-1,2))
    @test func(d0, 1.0) == 2.0
    # 
    d1 = fixpars(d0, (:a,))
    @test collectpars(d1) === ∅
    @test func(d1, 1.0) == 2.0
    # 
    d1′ = fixpars(d0, [:a])
    @test collectpars(d1′) === ∅
    @test func(d1′, 1.0) == 2.0
    # 
    d2 = fixpars(d0, (a=2,))
    @test collectpars(d2) === ∅
    @test func(d2, 1.0) == 3.0
end

@testset "example of usage" begin
    # tests
    snl = pdf(@. (x;p) -> exp(-(x-p.μ)^2/(2*p.σ^2)); p = (μ=1.4, σ=0.15), lims=(0, 3))
    bkg = pdf(@. (x;p) -> sqrt(x)*exp(-p.α*x); p = (α=1.3,), lims=(0, 3))
    bkg *= (fb=2.5,)
    pdf_sum = snl + bkg
    
    # generating
    data = generate(10000, pdf_sum; p=collectpars(pdf_sum));
    @test length(data) == 10000
    
    # fitting
    fr = fit_llh(data, pdf_sum)
    pfr = v2p(minimizer(fr), pdf_sum)
    @test length(pfr) == 4
    # 
    @test minimum(fr) ≤ llh(data, pdf_sum)
    # 
    invH_bfgs = invH(fr)
    @test size(invH_bfgs) == (4,4)
    @test cov(fr) == invH_bfgs
    @test length(errors(fr)) == 4
    # 
    invH_fd = invexacthessian(fr)
    rel_error_max = max((abs.(invH_bfgs .- invH_fd) ./ abs.(invH_fd))...)
    @test rel_error_max < 5e-2 # 5%
end

@testset "parameters to values" begin
    d = pdf(@. (e;p)->e^2+p.a; p=(a=1.0,), lims=(-1,2))
    @test p2v(d) == [1.0]
    @test p2v((a=3.0,), d) == [3.0]
end

@testset "fixed-shape pdf" begin
    d1 = fixedshapepdf(x->exp.(-(x .* 4).^2), (-1, 2))
    @test length(collectpars(d1)) == 0
end

g(x) = exp.(-(x .* 4).^2)
e(x) = exp.(-x)
mylims = (-1, 2)
# 
sum0 = sumpdf(g,e,mylims)
pdf1 = fixedshapepdf(g, mylims)
pdf2 = fixedshapepdf(e, mylims)
# 
sum1 = sumpdf(pdf1, pdf2)
sum2 = sumpdf(pdf1, pdf2, :xf)
sum3 = sumpdf(pdf1, pdf2, 0.3)
# 
xr = mylims[1]+rand()*(mylims[2]-mylims[1])

@testset "sum pdf" begin
    # 
    @test length(collectpars(sum0)) == 1
    @test length(collectpars(sum1)) == 1
    @test sum0(xr) ≈ sum1(xr)
    @test keys(collectpars(sum2))[1] == :xf
    @test length(collectpars(sum3)) == 0
end

@testset "example with sum pdf" begin
    ds = generate(5000, sum1)
    ft = fit_llh(ds,sum1; init_pars=[0.3])
    pfr = minimizer(ft)
    sum1_fit = fixpars(sum1, v2p(pfr,sum1))
    println("δf = ", (pfr[1] - sum1.p[1]) / sum1.p[1])
    @test (pfr[1] - sum1.p[1]) / sum1.p[1] < 0.05
end

@testset "no-parameters f and no-parameters normalized f" begin
    d0 = pdf(@. (e;p)->e^2+p.a; p=(a=1.0,), lims=(-1,2))

    f = noparsf(d0)
    xr = lims(d0)[1]+rand()*(lims(d0)[2]-lims(d0)[1])
    @test func(d0,xr) ≈ f(xr)
    #
    ananorm = ((8+1)/3+1.0*3)
    @test d0(xr) ≈ f(xr)/ananorm
end

sWeights_signal, sWeights_backgr = sWeights(pdf1, pdf2, 0.9)
xv = range(lims(pdf1)...,length=100)
sum_of_w = sWeights_signal(xv) + sWeights_backgr(xv)

@testset "sWeights" begin
    xv = range(lims(pdf1)...,length=100)
    sum_of_w = sWeights_signal(xv) + sWeights_backgr(xv)
    @test prod(sum_of_w .- sum_of_w[30] .< 1e-10)
end

@testset "cross-product PDF" begin
    # test
    pdf1 = pdf((x;p)->x.^2; lims=(-1,2), p=∅)
    pdf2 = pdf((x;p)->x.^4; lims=(-1,2), p=∅)
    X = xProductPDF(x=pdf1, y=pdf2)
    s = generate(100, X)
    # 
    @test length(s) == 100
    @test hasproperty(s[1], :x)
    @test hasproperty(s[1], :y)
    # 
    s = generate(50, X; Nbins=300)
    @test length(s) == 50
end

g1 = pdf(@. (x;p)->1/p.σ1*exp(-(x-p.μ1)^2/(2*p.σ1^2)); p=(μ1= 2.1, σ1=0.7 ), lims=(-3, 3))
g2 = pdf(@. (x;p)->1/p.σ2*exp(-(x-p.μ2)^2/(2*p.σ2^2)); p=(μ2=-0.7, σ2=0.7 ), lims=(-3, 3))
mm0 = MixedModel([g1, g2], (f1=0.33,))

@testset "parameters of the mixed model" begin
    @test length(collectpars(fixpars(mm0, (σ1=1.1,)))) == 4
    @test length(collectpars(fixpars(mm0, (:σ1,)))) == 4
    @test length(collectpars(fixpar(mm0, :σ1))) == 4
    @test length(collectpars(fixpar(mm0, :σ1, 1.1))) == 4
    #
end

@testset "Mixed Model" begin

    @test typeof(collectpars(mm0)) <: NamedTuple
    @test npars(mm0) == 5
    @test v2p(p2v(mm0),mm0) == collectpars(mm0)
    #
    sample = vcat((0.5.*randn(1000) .- 1.0), (0.7.*randn(100) .+ 2.0))
    sample = filter(x->-3<x<3, sample)
    # 
    fr = fit_llh(sample, mm0, init_pars=p2v(mm0))
    pfr = v2p(minimizer(fr), mm0)
    @test abs(pfr.μ2 + 1.0) < 0.1 && abs(pfr.μ1 - 2.0) < 0.5
    #
end

@testset "ensities" begin
    pdf1 = aGauss((mΩb = 6030, σ=17.0), (5600, 6400))
    @test collectpars(pdf1) === (mΩb = 6030, σ=17.0)
    pdf2 = aBreitWigner((mΩb = 6030, Γ=17.0), (5600, 6400))
    @test collectpars(pdf2) === (mΩb = 6030, Γ=17.0)
    pdf3 = aExp((τ = -1.1,), (-2, 2))
    @test collectpars(pdf3) === (τ = -1.1,)
    # 
    pdf4 = aDoubleGaussFixedRatio((m = 0.77, Γ=0.15), (0, 1.0); fixpars=(r=0.8,n=3))
    @test pdf4(0.77) != 0.0
    pdf5 = aBreitWignerConvGauss((m = 0.77, Γ=0.15), (0, 1.0); fixpars=(σ=0.03,))
    @test pdf5(0.77) != 0.0
    # 
    xv = range(-π, 2π, length=40)
    yv = map(x->x*cos(3x) - 3*sin(x), xv)
    pdf6 = aTabulated(xv,yv,(-π,π))
    @test length(collectpars(pdf6)) == 0
    @test pdf6(1.1) != 0.0
    @test pdf6(3π) == 0.0
end

@testset "plotting" begin
    scaletobinneddata(10,(0,1),10) ≈ 1.0
    scaletobinneddata(10, range(0,1,length=11)) ≈ 1.0
end