
using AlgebraPDF
using Test

g1 = FGauss((μ1=1.1,σ1=0.1))
g2 = FGauss((μ2=2.2,σ2=0.2))
g3 = FGauss((μ3=3.3,σ3=0.3))
d = AlgebraPDF.SumFunc([g1,g2,g3], (α1=0.1,α2=0.2,α3=0.3))
d2 = AlgebraPDF.SumFunc([g1,g2,g3], Ext(α1=0.1,α2=0.2,α3=0.3)) |> x->fixpar(x, :α1)

@testset "Sum of regular functions" begin 
    @test length(d) == 3
    @test abs(func(d, 1.5)) < 0.01
    @test abs(func(d, 1.5; p=pars(d)+(μ1=1.5,))) > 0.01
    @test func(d, 1.5; p=pars(d)+(α1=1500,)) > 0.01
    # 
    @test abs(d(1.1)-d(2.2)) < 0.01
    @test abs(d(2.2)-d(3.3)) < 0.01

    keys(fixedpars(d2)) == (:α1,)
    @test d2(1.1) == d(1.1)
end
#

ng1 = Normalized(g1, (1,4))
ng2 = Normalized(g2, (1,4))
ng3 = Normalized(g3, (1,4))
nd = AlgebraPDF.SumFunc([ng1,ng2,ng3], (α1=0.1,α2=0.2,α3=0.3))
# 
@testset "Sum of normalized functions" begin 

    @test typeof(d) <: AlgebraPDF.SumFunc{T} where T<:AbstractFunctionWithParameters
    @test typeof(nd) <: AlgebraPDF.SumFunc{T} where T<:AbstractPDF

    @test func(nd, 1.1) isa Number
    @test func(nd, [1.1,1,1])  isa Vector
    @test func(nd, 1:0.1:2)  isa Vector

    @test lims(nd) == lims(ng1)

    @test nd(1.1) > d(1.1) # due to normalization
end