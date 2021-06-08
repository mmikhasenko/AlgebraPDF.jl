
@testset "Operations on NamedTuple" begin
    p = (μ = 1.1, σ = 2.2, f = 2.2)

    sp = p - (:μ, :σ)
    @test keys(sp) == (:f,)
    @test sp == -(p, [:μ, :σ])
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
    ps0 = AlgebraPDF.TwoNamedTuples((a=1.1,b=2.2,c=3.3))
    d0 = FunctionWithParameters(f=(x;p)->x^2,p=ps0)
    # #
    d1 = fixpar(d0, :a)
    @test length(fixedpars(d1)) == 1 && length(freepars(d1)) == 2
    d2 = fixpars(d0, (a=3.3, b=6.6))
    @test pars(d2).a == 3.3 && pars(d2).b == 6.6
    # 
    @test releasepar(fixpar(d0, :c), :c) == d0
    # #
    @test length(freepars(d0)) == 3
    @test length(freepars(d1)) == 2
    #
    @test pars(updatepars(d0, (a=5.5,))).a == 5.5
end

