using AlgebraPDF
using Test
using SpecialFunctions

@testset "Basic operations" begin
    BW(s, m, Γ) = 1 / (m^2 - s - 1im*m*Γ)
    #
    pdf1 = pdf(@. (e;p)->abs2(BW(e^2, p.m1, p.Γ1));
        p0 = (m1=0.25, Γ1=2e-3), lims = (0, 0.15))
    #
    @test pdf1(rand()) != 0.0
    @test length(pdf1(rand(10))) == 10
    #
    pdf2 = pdf(@. (e;p)->abs2(BW(e^2, p.m2, p.Γ2));
        p0 = (m2=0.1, Γ2=14e-3), lims = (0, 0.15))
    #
    x0 = 1.1; v0 = rand(2)
    @test pdf2(x0, v0) == pdf2(x0; p=v2p(v0, pdf2))
    #
    pdf2 *= (f2=3.0,)
    @test length(pdf2.p0) == 3
    #
    pdf_sum = pdf1 + pdf2
    @test npars(pdf_sum) == 5
end

function df(f1,f2,lims; Ns = 10)
    xv = range(lims..., length=Ns)
    dyv = f1.(xv) .- f2.(xv)
    return sqrt(sum(abs2, dyv) / Ns)
end

@testset "convolution with gauss" begin

    @test AlgebraPDF.g(0,1) ≈ 1/sqrt(2π)
    @test AlgebraPDF.g(2,2) ≈ exp(-1/2)/sqrt(2π*4)
    # tests
    σ0 = 0.3
    aconv(e) = (erf(e/sqrt(2*σ0^2))+1)/2
    #
    nconv(e) = conv_with_gauss(e, x->x>0, σ0)
    @test df(nconv, aconv, (-1, 1); Ns=1000) < 0.01
    #
    step = pdf(@. (e;p)->e>0; p0=NamedTuple(), lims=(-1,1))
    smeared_step = conv_with_gauss(step, σ0)
    nconv(e) = smeared_step.f(e; p=smeared_step.p0)
    @test df(nconv, aconv, (-1, 1); Ns=1000) < 0.01
    # 
    smeared_step_sampling = conv_with_gauss_sampling(step, σ0; Ns=50)
    nconv(e) = smeared_step_sampling.f(e; p=smeared_step_sampling.p0)
    @test df(nconv, aconv, (-1, 1); Ns=1000) < 0.01
    # 
    smeared_step_sampling = conv_with_gauss_sampling(step, σ0; Ns=10)
    nconv(e) = smeared_step_sampling.f(e; p=smeared_step_sampling.p0)
    @test 0.01 < df(nconv, aconv, (-1, 1); Ns=1000) < 0.04
end

@testset "fix parameters example" begin
    d0 = pdf(@. (e;p)->e^2+p.a; p0=(a=1.0,), lims=(-1,2))
    @test d0.f(1.0; p=d0.p0) == 2.0
    # 
    d1 = fix_parameters(d0, (:a,))
    @test d1.p0 === NamedTuple()
    @test d1.f(1.0; p=d0.p0) == 2.0
    # 
    d1′ = fix_parameters(d0, [:a])
    @test d1′.p0 === NamedTuple()
    @test d1′.f(1.0; p=d0.p0) == 2.0
    # 
    d2 = fix_parameters(d0, (a=2,))
    @test d2.p0 === NamedTuple()
    @test d2.f(1.0; p=d2.p0) == 3.0
end

@testset "example of usage" begin
    # tests
    snl = pdf(@. (x;p) -> exp(-(x-p.μ)^2/(2*p.σ^2)); p0 = (μ=1.4, σ=0.15), lims=(0, 3))
    bkg = pdf(@. (x;p) -> sqrt(x)*exp(-p.α*x); p0 = (α=1.3,), lims=(0, 3))
    bkg *= (fb=2.5,)
    pdf_sum = snl + bkg
    
    # generating
    data = generate(10000, pdf_sum; p=pdf_sum.p0);
    @test length(data) == 10000
    
    # fitting
    fr = fit_llh(data, pdf_sum; init_pars=p2v(pdf_sum.p0, pdf_sum))
    pfr = v2p(fr.minimizer, pdf_sum)
    @test length(pfr) == 4
end

@testset "parameters to values" begin
    d = pdf(@. (e;p)->e^2+p.a; p0=(a=1.0,), lims=(-1,2))
    @test p2v(d) == [1.0]
    @test p2v((a=3.0,), d) == [3.0]
end