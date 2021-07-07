using AlgebraPDF
using Test

@testset "Fit gauss with Minuit" begin 
    d0 = Normalized(FGauss((μ=1.1, σ=0.3)), (-3, 5))
    data = randn(1000)
    fs0 = fit(d0, data, MigradAndHesse())

    @test keys(fs0.parameters) == (:μ,:σ)
    @test typeof(Tuple(fs0.measurements)) <: Tuple{Vararg{Measurement}}

    d1 = fixpar(Normalized(FGauss(Ext(μ=1.1, σ=0.9)), (-3, 5)), :σ)
    fs1 = fit(d1, data, MigradAndHesse())

    @test keys(fs1.parameters) == (:μ,)
    @test pars(fs1.best_model).σ == 0.9
end


@testset "Fit gauss with Optim" begin 
    d0 = Normalized(FGauss((μ=0.01, σ=0.9)), (-3, 5))
    data = randn(1000)
    fs0 = fit(d0, data, OptimHesseApprox())

    @test keys(fs0.parameters) == (:μ,:σ)
    @test typeof(Tuple(fs0.measurements)) <: Tuple{Vararg{Measurement}}

    d1 = fixpar(Normalized(FGauss(Ext(μ=0.01, σ=0.9)), (-3, 5)), :σ)
    fs1 = fit(d1, data, OptimHesseApprox())

    @test keys(fs1.parameters) == (:μ,)
    @test pars(fs1.best_model).σ == 0.9
end