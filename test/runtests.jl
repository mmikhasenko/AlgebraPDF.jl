using AlgebraPDF
using StaticArrays
using Test
using SpecialFunctions
using LinearAlgebra
using Measurements


@testset "Operations on NamedTuple" begin
    p = (μ = 1.1, σ = 2.2, f = 2.2)

    sp = p - (:μ, :σ)
    @test keys(sp) == (:f,)
    @test sp == -(p, [:μ, :σ])
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
    # 
    @test (a=1,b=2) + (d=1,c=2) == (a=1,b=2,d=1,c=2)
    @test (a=1,b=2,d=1,c=2) - (d=1,) == (a=1,b=2,c=2)
    @test (a=1,b=2,d=1,c=2) - :d == (a=1,b=2,c=2)
    @test (a=1,b=2,d=1,c=2) - [:d, :a] == (b=2,c=2)
    @test (a=1,b=2,d=1,c=2) - (:d, :a) == (b=2,c=2)
end

@testset "parameter logic" begin
    # 
    @test nt(:d) == (d=0.0, )
    @test nt(:d, 3) == (d=3, )
    @test nt(:d, (1.1,0.1)) == (d = (1.1,0.1), )
    #
    ps0 = AlgebraPDF.Parameters((a=1.1,b=2.2,c=3.3))
    #
    ps1 = fixpar(ps0, :a)
    @test length(fixedpars(ps1)) == 1 && length(freepars(ps1)) == 2
    ps2 = fixpars(ps0, (a=3.3, b=6.6))
    @test ps2.a == 3.3 && ps2.b == 6.6
    # 
    @test releasepar(fixpar(ps0, :c), :c) == ps0
    # 
    ps3 = constrainpar(ps0, :a, 1.1, 0.1)
    @test length(constrainedpars(ps3)) == 1 && length(freepars(ps3)) == 3
    @test unconstrainpar(ps3, :a) == ps0
    #
    @test length(freepars(ps0)) == 3
    @test length(freepars(ps1)) == 2
    @test length(freepars(ps3)) == 3
    #
    @test updatepars(ps0, (a=5.5,)).a == 5.5
end

@testset "Parameters of pdf" begin
    d = pdf(@. (e;p)->e^2+p.a; p=(a=1.0,), lims=(-1,2))
    @test typeof(pars(d)) <: AlgebraPDF.Parameters
    @test length(freepars(d)) == 1
    @test length(fixedpars(d)) == 0
end

@testset "fix parameters example" begin
    d0 = pdf(@. (e;p)->e^2+p.a; p=(a=1.0,), lims=(-1,2))
    @test d0(1.0) ≈ 2.0/(9/3+3*pars(d0).a)
    #
    d1 = fixpar(d0, :a, 2.9)
    @test d1(1.0) ≈ 3.9/(9/3+3*pars(d1).a)
    # 
    d1′ = fixpars(d0, (:a,))
    @test d1′ == fixpars(d0, [:a])
    @test d1′ == fixpars(d0, (a=1.0,))
end

@testset "parameters to values" begin
    d = pdf(@. (e;p)->e^2+p.a; p=(a=1.0,), lims=(-1,2))
    @test p2v(d) == [1.0]
    @test p2v((a=3.0,), d) == [3.0]
end

@testset "fixed-shape pdf" begin
    d1 = fixedshapepdf(x->exp.(-(x .* 4).^2), (-1, 2))
    @test length(freepars(d1)) == 0
end

@testset "no-parameters f and no-parameters normalized f" begin
    d0 = pdf(@. (e;p)->e^2+p.a; p=(a=1.0,), lims=(-1,2))

    f = noparsf(d0)
    xr = lims(d0)[1]+rand()*(lims(d0)[2]-lims(d0)[1])
    @test func(d0,xr;p=freepars(d0)) ≈ f(xr)
    #
    ananorm = ((8+1)/3+1.0*3)
    @test d0(xr) ≈ f(xr)/ananorm
end


@testset "example of usage" begin
    # tests
    snl = pdf(@. (x;p) -> exp(-(x-p.μ)^2/(2*p.σ^2)); p = (μ=1.4, σ=0.15), lims=(0, 3))
    bkg = pdf(@. (x;p) -> sqrt(x)*exp(-p.α*x); p = (α=1.3,), lims=(0, 3))
    bkg *= (fb=2.5,)
    pdf_sum = snl + bkg
    
    # generating
    data = generate(10000, pdf_sum; p=freepars(pdf_sum));
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
    @test covmat(fr) == invH_bfgs
    @test length(errors(fr)) == 4
    # 
    invH_fd = invexacthessian(fr)
    rel_error_max = max((abs.(diag(invH_bfgs) .- diag(invH_fd)) ./ abs.(diag(invH_fd)))...)
    @test rel_error_max < 10e-2 # 10%
    # 
    mfr = measurements(fr)
    @test typeof([measurements(fr)...]) <: Vector{Measurement{T}} where T
    mfr_exact = measurements(fr, true)
    @test Measurements.uncertainty.(mfr) == errors(fr)
    @test Measurements.uncertainty.(mfr_exact) == sqrt.(diag(invH_fd))
end


@testset "pdf with NamedTuple" begin
    d = pdf((x; p)->x.*5 .+ p.a, (-2, 3), 1, (a=1.1, ))
    # 
    @test d(1.1) != 0.0
    # 
    @test_throws DomainError fixpar(d, :a, 1.2)
    @test_throws DomainError fixpars(d, (a=1.2,))
    @test_throws DomainError releasepar(d, :a)
    @test_throws DomainError constrainpar(d, :a, 1.1, 0.2)
    @test_throws DomainError unconstrainpar(d, :a)
    # 
end



#              _|_|_|                                
#    _|_|_|  _|          _|_|_|  _|    _|    _|_|_|  
#  _|    _|  _|  _|_|  _|    _|  _|    _|  _|_|      
#  _|    _|  _|    _|  _|    _|  _|    _|      _|_|  
#    _|_|_|    _|_|_|    _|_|_|    _|_|_|  _|_|_|    


@testset "Precodded pdfs" begin
    pdf1 = aGauss((mΩb = 6030, σ=17.0), (5600, 6400))
    @test freepars(pdf1) === (mΩb = 6030, σ=17.0)
    pdf2 = aBreitWigner((mΩb = 6030, Γ=17.0), (5600, 6400))
    @test freepars(pdf2) === (mΩb = 6030, Γ=17.0)
    pdf3 = aExp((τ = -1.1,), (-2, 2))
    @test freepars(pdf3) === (τ = -1.1,)

    pdfPowExp = aPowExp((n=3.3, τ=-4.0), (0, 2.2))
    @test pdfPowExp(1.1) != 0.0
    @test freepars(pdfPowExp) === (n=3.3, τ=-4.0)
    # 
    pars3 = (c0=1.1, c1=2.2, c2=3.3, c4=4.4)
    pdfPpol3 = aPol(pars3, (-2,12))
    @test func(pdfPpol3, 1.1; p=pars3) == sum(1.1^(i-1)*c for (i,c) in enumerate(pars3))
    @test freepars(pdfPpol3) === pars3
    # 
    pdf4 = aDoubleGaussFixedRatio((m = 0.77, Γ=0.15), (0, 1.0); fixpars=(r=0.8,n=3))
    @test pdf4(0.77) != 0.0
    pdf5 = aBreitWignerConvGauss((m = 0.77, Γ=0.15), (0, 1.0); fixpars=(σ=0.03,))
    @test pdf5(0.77) != 0.0
    # 
    xv = range(-π, 2π, length=40)
    yv = map(x->x*cos(3x) - 3*sin(x), xv)
    pdf6 = aTabulated(xv,yv,(-π,π))
    @test length(freepars(pdf6)) == 0
    @test pdf6(1.1) != 0.0
    @test pdf6(3π) == 0.0
    @test pdf6(-π) != 0.0
    @test pdf6(2π) != 0.0
end

# here is MWE of the normalization problem
# let 
#     d = aGauss((a=0.01,b = 1.1), (-3,3))
#     d(1.1, [36.1,0.1])
# end


#    _|_|_|    _|_|    _|_|_|    _|      _|  
#  _|        _|    _|  _|    _|  _|      _|  
#  _|        _|    _|  _|    _|    _|  _|    
#    _|_|_|    _|_|    _|    _|      _|      


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
    func(step, 1.1; p=∅)
    smeared_step = conv_with_gauss(step, σ0)
    nconv(e) = func(smeared_step, e; p=∅)
    # smeared_step(1.1)
    func(smeared_step, 1.1; p=∅)
    @test df(nconv, aconv, (-1, 1); Ns=1000) < 0.01
    # 
    smeared_step_sampling = conv_with_gauss_sampling(step, σ0; Ns=50)
    nconv(e) = func(smeared_step_sampling, e; p = freepars(smeared_step_sampling))
    @test df(nconv, aconv, (-1, 1); Ns=1000) < 0.01
    # 
    smeared_step_sampling = conv_with_gauss_sampling(step, σ0; Ns=10)
    nconv(e) = func(smeared_step_sampling, e; p = freepars(smeared_step_sampling))
    @test 0.01 < df(nconv, aconv, (-1, 1); Ns=1000) < 0.04
end

#                            _|    _|      _|                    _|  
#  _|_|_|  _|_|    _|    _|  _|  _|_|_|_|                    _|_|_|  
#  _|    _|    _|  _|    _|  _|    _|      _|  _|_|_|_|_|  _|    _|  
#  _|    _|    _|  _|    _|  _|    _|      _|              _|    _|  
#  _|    _|    _|    _|_|_|  _|      _|_|  _|                _|_|_|  


@testset "cross-product PDF" begin
    # test
    pdf1 = pdf((x;p)->x.^2; lims=(-1,2), p=∅)
    pdf2 = pdf((x;p)->x.^4; lims=(-1,2), p=∅)
    X = xProductPDF(;x=pdf1, y=pdf2)
    s = generate(100, X)
    # 
    @test length(s) == 100
    @test hasproperty(s[1], :x)
    @test hasproperty(s[1], :y)
    # 
    s = generate(50, X; Nbins=300)
    @test length(s) == 50
end


#                  _|                            _|                                  _|  
#  _|_|_|  _|_|        _|    _|    _|_|      _|_|_|  _|_|_|  _|_|      _|_|      _|_|_|  
#  _|    _|    _|  _|    _|_|    _|_|_|_|  _|    _|  _|    _|    _|  _|    _|  _|    _|  
#  _|    _|    _|  _|  _|    _|  _|        _|    _|  _|    _|    _|  _|    _|  _|    _|  
#  _|    _|    _|  _|  _|    _|    _|_|_|    _|_|_|  _|    _|    _|    _|_|      _|_|_|  


g1 = pdf(@. (x;p)->1/p.σ1*exp(-(x-p.μ1)^2/(2*p.σ1^2)); p=(μ1= 2.1, σ1=0.7 ), lims=(-3, 3))
g2 = pdf(@. (x;p)->1/p.σ2*exp(-(x-p.μ2)^2/(2*p.σ2^2)); p=(μ2=-0.7, σ2=0.7 ), lims=(-3, 3))
mm0 = MixedModel([g1, g2], (f1=0.33,))

@testset "Parameters of the mixed model" begin
    @test typeof(pars(mm0)) <: AlgebraPDF.Parameters
    @test length(freepars(mm0)) == 5
    @test length(fixedpars(mm0)) == 0
    # 
    @test fixpar(mm0, :σ1, 1.1) == fixpars(mm0, (σ1=1.1,))
    @test fixpars(mm0, (:σ1,)) == fixpar(mm0, :σ1) 
    mm_fixed = fixpar(mm0, :σ1, 1.1)
    @test length(freepars(mm_fixed)) == 4
    #
    mm′ = fixpar(mm0, :μ1, 2.2)
    @test length(fixedpars(fractions(mm′))) == 0
    @test length(fixedpars(mm′.components[1])) == 1
end

@testset "Integrals of the MM" begin
    @test integral(mm0, (-1,1)) < 1.0
    @test prod(integrals(mm0, (-1,1)) .< 1.0)
    #
    mm′ = fixpar(mm0, :μ1, 2.2)
    @test integral(mm′, (-1,1)) < 1.0
    @test prod(integrals(mm′, (-1,1)) .< 1.0)
end

@testset "Mixed Model" begin
    @test typeof(freepars(mm0)) <: NamedTuple
    @test npars(mm0) == 5
    @test v2p(p2v(mm0),mm0) == freepars(mm0)
    #
    sample = vcat((0.5 .* randn(1000) .- 1.0), (0.7.*randn(100) .+ 2.0))
    sample = filter(x->-3<x<3, sample)
    # 
    fr = fit_llh(sample, mm0, init_pars=p2v(mm0))
    pfr = v2p(minimizer(fr), mm0)
    @test abs(pfr.μ2 + 1.0) < 0.5
    @test abs(pfr.μ1 - 2.0) < 1.0
end

@testset "Generate mixed model" begin
    p1 = fixedshapepdf(x->(x .> 0), (-1,1))
    p2 = fixedshapepdf(x->0.01.*(x .< 0), (-1,1))
    mm = MixedModel([p1,p2],(f1 = 0.01,))
    Nd = 1000
    data = generate(Nd, mm)
    @test sum(x->x>0, data) < 100
end

@testset "Update parameters of the mixed model" begin
    g1 = aGauss((μ1=1.1, σ=0.3), (-2, 7))
    g2 = aGauss((μ2=2.1, σ=0.3), (-2, 7))
    # 
    mm = MixedModel([g1, g2], (f1=0.1,))
    mm = fixpar(mm, :f1, 0.3)
    mm = updatepars(mm, (f1=0.7,))
    @test fractions(mm).f1 == 0.7
end


#            _|              _|      _|      _|                      
#  _|_|_|    _|    _|_|    _|_|_|_|_|_|_|_|      _|_|_|      _|_|_|  
#  _|    _|  _|  _|    _|    _|      _|      _|  _|    _|  _|    _|  
#  _|    _|  _|  _|    _|    _|      _|      _|  _|    _|  _|    _|  
#  _|_|_|    _|    _|_|        _|_|    _|_|  _|  _|    _|    _|_|_|  
#  _|                                                            _|  
#  _|                                                        _|_|    


@testset "Plotting utils" begin
    scaletobinneddata(10,(0,1),10) ≈ 1.0
    scaletobinneddata(10, range(0,1,length=11)) ≈ 1.0
end

#                                            _|                
#    _|_|_|    _|_|    _|_|_|      _|_|_|  _|_|_|_|  _|  _|_|  
#  _|        _|    _|  _|    _|  _|_|        _|      _|_|      
#  _|        _|    _|  _|    _|      _|_|    _|      _|        
#    _|_|_|    _|_|    _|    _|  _|_|_|        _|_|  _|        


@testset "Constrained fit" begin
    @test AlgebraPDF.chi2((a=(1,0.2), b=(3,0.2)); p = (a=1.1, b=2.2)) ≈
        (0.1/0.2)^2+(0.8/0.2)^2
    #
    d = aGauss((a=0.01,b = 1.1), (-3,3))
    data = randn(1000)
    my_fr = fit_llh(data, d)
    my_pfr = v2p(minimizer(my_fr), d)
    # 
    constraints = (a = (0.2,0.01),)
    my_fr2 = fit_llh_with_constraints(data, d, constraints)
    my_pfr2 = v2p(minimizer(my_fr2), d)
    #
    @test abs(my_pfr2.a - constraints.a[1]) < abs(my_pfr.a - constraints.a[1])
end

#            _|                      
#  _|_|_|    _|  _|    _|    _|_|_|  
#  _|    _|  _|  _|    _|  _|_|      
#  _|    _|  _|  _|    _|      _|_|  
#  _|_|_|    _|    _|_|_|  _|_|_|    
#  _|                                
#  _|                                


@testset "Arithmetic operations" begin
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
    @test length(freepars(pdf2)) == 3
    #
    pdf_sum = pdf1 + pdf2
    @test npars(pdf_sum) == 5

    pdf_ratio = pdf1 / pdf2
    @test npars(pdf_ratio) == 5
end

                                                                                   
#                                                      _|      _|_|  
#    _|_|_|  _|    _|  _|_|_|  _|_|    _|_|_|      _|_|_|    _|      
#  _|_|      _|    _|  _|    _|    _|  _|    _|  _|    _|  _|_|_|_|  
#      _|_|  _|    _|  _|    _|    _|  _|    _|  _|    _|    _|      
#  _|_|_|      _|_|_|  _|    _|    _|  _|_|_|      _|_|_|    _|      
#                                      _|                            
#                                      _|                            

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

@testset "Sum pdf" begin
    @test length(freepars(sum0)) == 1
    @test length(freepars(sum1)) == 1
    @test sum0(xr) ≈ sum1(xr)
    @test keys(freepars(sum2))[1] == :xf
    @test length(freepars(sum3)) == 0
end

@testset "Example with sum pdf" begin
    ds = generate(5000, sum1)
    ft = fit_llh(ds,sum1; init_pars=[0.3])
    pfr = minimizer(ft)
    sum1_fit = fixpars(sum1, v2p(pfr,sum1))
    println("δf = ", (pfr[1] - freepars(sum1)[1]) / freepars(sum1)[1])
    @test (pfr[1] - freepars(sum1)[1]) / freepars(sum1)[1] < 0.05
end


#            _|_|_|    _|              _|      
#    _|_|_|  _|    _|  _|    _|_|    _|_|_|_|  
#  _|_|      _|_|_|    _|  _|    _|    _|      
#      _|_|  _|        _|  _|    _|    _|      
#  _|_|_|    _|        _|    _|_|        _|_|  


sWeights_signal, sWeights_backgr = sWeights(pdf1, pdf2, 0.9)
xv = range(lims(pdf1)...,length=100)
sum_of_w = sWeights_signal(xv) + sWeights_backgr(xv)

@testset "sWeights" begin
    xv = range(lims(pdf1)...,length=100)
    sum_of_w = sWeights_signal(xv) + sWeights_backgr(xv)
    @test prod(sum_of_w .- sum_of_w[30] .< 1e-10)
end
