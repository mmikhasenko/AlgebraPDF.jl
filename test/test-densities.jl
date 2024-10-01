using AlgebraPDF
using Test
using AlgebraPDF.QuadGK

#              _|_|_|
#    _|_|_|  _|          _|_|_|  _|    _|    _|_|_|
#  _|    _|  _|  _|_|  _|    _|  _|    _|  _|_|
#  _|    _|  _|    _|  _|    _|  _|    _|      _|_|
#    _|_|_|    _|_|_|    _|_|_|    _|_|_|  _|_|_|



@testset "Standard functions" begin
    @test AlgebraPDF.standardgauss(0,1) ≈ 1/sqrt(2π)
    @test AlgebraPDF.standardgauss(2,2) ≈ exp(-1/2)/sqrt(2π*4)
    x, m, Γ = 1.0, 1.0, 0.1
    @test AlgebraPDF.amplitudeBW(x, m, Γ) ≈ 1im
    @test AlgebraPDF.amplitudeBWsq(x, m, Γ) ≈ 1.0
    #
    x, σ, r, n = 0, 0.1, 0.3, 5
    @test AlgebraPDF.standarddoublegauss(x,σ,r,n) ≈ 1/sqrt(2π)/σ*(r+(1-r)/n)
    #
end


@testset "Precodded pdfs" begin
    dGauss = Normalized(FGauss((μ1=1.0,σ1=0.3)), (-2,2))
    @test keys(pars(dGauss)) == (:μ1,:σ1)
    exp_gauss_at_zero = 1/sqrt(2π)/pars(dGauss).σ1
    @test abs(dGauss(1.0) - exp_gauss_at_zero) / exp_gauss_at_zero < 1e-3
end

@testset "Precodded pdfs" begin
    d2 = FBreitWigner((mΩb = 6030, Γ=17.0))
    @test freepars(d2) == (mΩb = 6030, Γ=17.0)
    @test real(d2(freepars(d2).mΩb)) ≈ 0
    #
    d3 = FExp((τ = -1.1,))
    @test freepars(d3) === (τ = -1.1,)
    #
    d4 = FPowExp((n=3.3, τ=-4.0))
    @test d4(1.1) != 0.0
    @test freepars(d4) == (n=3.3, τ=-4.0)
    #
    pars3 = (c0=1.1, c1=2.2, c2=3.3, c4=4.4)
    d5 = FPol(pars3)
    @test func(d5, 1.1; p=pars3) == sum(1.1^(i-1)*c for (i,c) in enumerate(pars3))
    @test freepars(d5) == pars3
    #
    d6 = FDoubleGaussFixedRatio((m = 0.77, Γ=0.15, r=0.8, n=3))
    @test d6(0.77) != 0.0
    #
    d7 = FBreitWignerConvGauss((m = 0.77, Γ=0.15, σ=0.03))
    @test d7(0.77) != 0.0
    #
    xv = range(-π, 2π, length=40)
    yv = map(x->x*cos(3x) - 3*sin(x), xv)
    d8 = FTabulated(xv,yv)
    @test length(freepars(d8)) == 0
    @test d8(1.1) != 0.0
    @test d8( 3π) == 0.0
    @test d8( -π) != 0.0
    @test d8( 2π) != 0.0
    #
end

# @testset "Crystal Ball" begin
let
    cb1 = (μ=0.0, σ=1, α=2.0, n=2) |> FLeftSideCrystalBall
    a, b = cb1.(-cb1.p.α .+ 1e-5 .* [-1,1])
    # match on the connection point
	@test isapprox(a,b; rtol = 1e-4)

    cb2 = (μ=0.0, σ=1, α=2.0, n=2) |> FRightSideCrystalBall
    a, b = cb2.(cb2.p.α .+ 1e-5 .* [1,-1])
    # match on the connection point
	@test isapprox(a,b; rtol = 1e-4)

    cb3 = (μ=0.0, σ=1, α=2.0, n=2) |> FDoubleSideCrystalBall
    a, b = cb3.(-cb3.p.α .+ 1e-5 .* [-1,1])
    # match on the connection point
	@test isapprox(a,b; rtol = 1e-4)
    # symmetric
    @test  prod((a, b) .≈ cb3.(cb3.p.α .+ 1e-5 .* [1,-1]))
    #
    @test quadgk(cb1, -Inf, Inf)[1] ≈ 1.0
    @test quadgk(cb2, -Inf, Inf)[1] ≈ 1.0
    @test quadgk(cb3, -Inf, Inf)[1] ≈ 1.0
end
