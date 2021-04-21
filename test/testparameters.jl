
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
