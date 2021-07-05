
using AlgebraPDF

@testset "Parameter utils" begin
    @test nt(:d) == (d=0.0, )
    @test nt(:d, 3) == (d=3, )
    @test nt(:d, (1.1,0.1)) == (d = (1.1,0.1), )
end

@testset "Operations on NamedTuple" begin
    @test (a=1,b=2) + (d=1,c=2) == (a=1,b=2,d=1,c=2)
    @test (a=1,b=2,d=1,c=2) - (d=1,) == (a=1,b=2,c=2)
    @test (a=1,b=2,d=1,c=2) - :d == (a=1,b=2,c=2)
    @test (a=1,b=2,d=1,c=2) - [:d, :a] == (b=2,c=2)
    @test (a=1,b=2,d=1,c=2) - (:d, :a) == (b=2,c=2)
end


@testset "Parameter structure: NamedTuple" begin
    p0 = (a=1,b=2)
    @test pars(p0) == p0
    @test freepars(p0) == p0
    @test fixedpars(p0) == ∅
    @test AlgebraPDF.updatevalueorflag(p0, :a, true, 2) == (a=2,b=2)
    @test_throws ArgumentError AlgebraPDF.updatevalueorflag(p0, :a, false)
end


@testset "Parameter structure: NamedTuple" begin
    #
    p0 = FlaggedNamedTuple(a=1,b=2)
    @test p0.a == 1
    @test p0.b == 2
    # 
    p1 = AlgebraPDF.updatevalueorflag(p0, :a, true, 2)
    @test pars(p1) == (a=2,b=2)
    #
    p2 = AlgebraPDF.updatevalueorflag(p0, :a, false)
    @test pars(p0) == pars(p2)
    @test freepars(p2) == (b=2,)
    @test fixedpars(p2) == (a=1,)
    #
    @test p0 == FlaggedNamedTuple((a=1,b=2))
    @test p2 == FlaggedNamedTuple(p2)
    # 
    @test keys(p1) == keys(p2)
    @test keys(p2, true) == (:b,)
    @test keys(p2, false) == (:a,)
end
