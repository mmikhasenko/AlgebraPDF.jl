

#                  _|                            _|                                  _|  
#  _|_|_|  _|_|        _|    _|    _|_|      _|_|_|  _|_|_|  _|_|      _|_|      _|_|_|  
#  _|    _|    _|  _|    _|_|    _|_|_|_|  _|    _|  _|    _|    _|  _|    _|  _|    _|  
#  _|    _|    _|  _|  _|    _|  _|        _|    _|  _|    _|    _|  _|    _|  _|    _|  
#  _|    _|    _|  _|  _|    _|    _|_|_|    _|_|_|  _|    _|    _|    _|_|      _|_|_|  


g1 = pdf((x;p)->1/p.σ1*exp(-(x-p.μ1)^2/(2*p.σ1^2)); p=TwoNamedTuples((μ1= 2.1, σ1=0.7 )), lims=(-3, 3))
g2 = pdf((x;p)->1/p.σ2*exp(-(x-p.μ2)^2/(2*p.σ2^2)); p=TwoNamedTuples((μ2=-0.7, σ2=0.7 )), lims=(-3, 3))
mm0 = MixedModel([g1, g2], (f1=0.33,))
mm1 = MixedModel([g1, g2], AlgebraPDF.TwoNamedTuples((f1=0.33,)))

@testset "Parameters of the mixed model" begin
    @test length(freepars(mm1)) == 5
    @test length(fixedpars(mm1)) == 0
    # 
    @test_throws ArgumentError fixpar(mm0, :f1, 1.1)
    @test fixpar(mm1, :σ1, 1.1) == fixpars(mm1, (σ1=1.1,))
    @test fixpars(mm1, (:σ1,)) == fixpar(mm1, :σ1) 
    mm_fixed = fixpar(mm1, :σ1, 1.1)
    @test length(freepars(mm_fixed)) == 4
    #
    mm′ = fixpar(mm1, :μ1, 2.2)
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

# @testset "Mixed Model" begin
#     @test typeof(freepars(mm0)) <: NamedTuple
#     @test nfreepars(mm0) == 5
#     @test v2p(p2v(mm0),mm0) == freepars(mm0)
#     #
#     sample = vcat((0.5 .* randn(1000) .- 1.0), (0.7.*randn(100) .+ 2.0))
#     sample = filter(x->-3<x<3, sample)
#     # 
#     fr = fit_llh(sample, mm0, init_pars=p2v(mm0))
#     pfr = v2p(minimizer(fr), mm0)
#     @test abs(pfr.μ2 + 1.0) < 0.5
#     @test abs(pfr.μ1 - 2.0) < 1.0
# end

# @testset "Generate mixed model" begin
#     p1 = fixedshapepdf(x->(x>0), (-1,1))
#     p2 = fixedshapepdf(x->0.01*(x<0), (-1,1))
#     mm = MixedModel([p1,p2],(f1 = 0.01,))
#     Nd = 1000
#     data = generate(Nd, mm)
#     @test sum(x->x>0, data) < 100
# end

@testset "Update parameters of the mixed model" begin
    g2 = PDFWithParameters(FGauss(Ext(μ2=2.1, σ=0.3)), (-2, 7))
    g1 = PDFWithParameters(FGauss(Ext(μ1=1.1, σ=0.3)), (-2, 7))
    # 
    mm = MixedModel([g1, g2], TwoNamedTuples((f1=0.1,)))
    mm = fixpar(mm, :f1, 0.3)
    mm = updatepars(mm, (f1=0.7,))
    @test fractions(mm).f1 == 0.7
end
