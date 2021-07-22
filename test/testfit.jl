using AlgebraPDF
using Test
using Random
using Measurements

@testset "Fit gauss with Minuit" begin 
    d0 = Normalized(FGauss((μ=1.1, σ=0.3)), (-3, 5))
    Random.seed!(1234)
    data = randn(1000)
    fs0 = fit(d0, data, MigradAndHesse())

    @test keys(fs0.parameters) == (:μ,:σ)
    @test typeof(Tuple(fs0.measurements)) <: Tuple{Vararg{Measurement}}

    d1 = fixpar(Normalized(FGauss(Ext(μ=1.1, σ=0.9)), (-3, 5)), :σ)
    fs1 = fit(d1, data, MigradAndHesse())

    @test keys(fs1.parameters) == (:μ,)
    @test pars(fs1.best_model).σ == 0.9
end

@testset "Extended Likelihood fit" begin 

    d = Normalized(FGauss((μ=1.1,σ=0.1)), (-2, 3))
    Random.seed!(4321)
    data = filter(x->inrange(x,lims(d)), randn(1000) .+ 0.5)

    fr = fit(d, data)
    @test 0.3 < fr.parameters.μ < 0.7
    @test 0.7 < abs(fr.parameters.σ) < 1.3

    s = FSum([d],(α=2.2,))

    @test integral(s) == 2.2
    @test s(1.1) == 2.2*d(1.1)
    @test s([1,2,3]) == 2.2*d([1,2,3])

    enll = Extended(NegativeLogLikelihood(s, data))
    fr2 = fit(enll)
    @test 0.3 < fr2.parameters.μ < 0.7
    @test 0.7 < abs(fr2.parameters.σ) < 1.3
    @test 900 < fr2.parameters.α < 1100
end
#
# Optim often fails with `isfinite(phi_c) && isfinite(dphi_c)` 
# removed from the tests
# 
@testset "Fit gauss with Optim" begin 
    d0 = Normalized(FGauss((μ=0.01, σ=0.9)), (-3, 5))
    Random.seed!(1234)
    data = randn(1000)
    fs0 = fit(d0, data, BFGSApproxHesse())

    @test keys(fs0.parameters) == (:μ,:σ)
    @test typeof(Tuple(fs0.measurements)) <: Tuple{Vararg{Measurement}}

    d1 = fixpar(Normalized(FGauss(Ext(μ=0.01, σ=0.9)), (-3, 5)), :σ)
    fs1 = fit(d1, data, BFGSApproxHesse())

    @test keys(fs1.parameters) == (:μ,)
    @test pars(fs1.best_model).σ == 0.9
end