using AlgebraPDF
using StaticArrays
using Test
using SpecialFunctions
using LinearAlgebra
using Measurements


@testset "func on scalars" begin
    @test func(1.1, 3.3) == 1.1
end

include("testparameters.jl")
include("testpdf.jl")
include("testdensities.jl")
include("testconvolution.jl")
include("testmultid.jl")
include("testmixedmodel.jl")

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
    pdf1 = pdf((e;p)->abs2(BW(e^2, p.m1, p.Γ1));
        p = (m1=0.25, Γ1=2e-3), lims = (0, 0.15))
    #
    @test pdf1(rand()) != 0.0
    @test length(pdf1(rand(10))) == 10
    #
    pdf2 = pdf((e;p)->abs2(BW(e^2, p.m2, p.Γ2));
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

g(x) = exp(-(4x)^2)
e(x) = exp(-x)
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

# @newfunc GG1(x;p) = x^2+p.a*x^3
# @newfunc GG2(x;p) = x^1+p.b*x^2
# h1 = GG1(p=(a=0.5,))
# h2 = GG2(p=(b=0.5,))

# h12 = h1+h2
# @testset "sum of functions" begin
#     @test func(h12,1.1) == func(h1,1.1) + func(h2,1.1)
#     @test keys(freepars(h12)) == (:a,:b,:α)
# end

# h1sq = abs2(h1)
# @testset "abs2 of functions" begin
#     @test pars(h1sq) == pars(h1)
#     @test func(h1sq, 1.1) == func(h1, 1.1)^2
# end

# @testset "implementation of copy" begin
#     @test func(h1, 1.1) ≈ func(copy(h1, pars(h1)), 1.1)
#     @test func(h12, 1.1) ≈ func(copy(h12, pars(h12)), 1.1)
#     @test func(h1sq, 1.1) ≈ func(copy(h1sq, pars(h1sq)), 1.1)
# end
