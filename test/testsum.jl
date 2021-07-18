
using AlgebraPDF
using Test

g1 = FGauss((μ1=1.1,σ1=0.1))
g2 = FGauss((μ2=2.2,σ2=0.2))
g3 = FGauss((μ3=3.3,σ3=0.3))
d = AlgebraPDF.FSum([g1,g2,g3], (α1=0.1,α2=0.2,α3=0.3))
d2 = AlgebraPDF.FSum([g1,g2,g3], Ext(α1=0.1,α2=0.2,α3=0.3)) |> x->fixpar(x, :α1)

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
nd = FSum([ng1,ng2,ng3], (α1=0.1,α2=0.2,α3=0.3))
# 
@testset "Sum of normalized functions" begin 

    @test typeof(d) <: FSum{T} where T<:AbstractFunctionWithParameters
    @test typeof(nd) <: FSum{T} where T<:AbstractPDF

    @test func(nd, 1.1) isa Number
    @test func(nd, [1.1,1,1])  isa Vector
    @test func(nd, 1:0.1:2)  isa Vector

    @test lims(nd) == lims(ng1)

    @test nd(1.1) > d(1.1) # due to normalization
end

@testset "Addition: FSum + FSum" begin
    s12 = FSum([g1, g2], (α1=1.1, α2=2.2))
    s3 = FSum([g3], (α3=3.3,))
    s123 = FSum([g1, g2, g3], (α1=1.1, α2=2.2, α3=3.3))

    @test s12 + s3 == s123

    s12_e = FSum([g1, g2], Ext(α1=1.1, α2=2.2))
    s3_e = FSum([g3], Ext(α3=3.3,))
    s123_e = FSum([g1, g2, g3], Ext(α1=1.1, α2=2.2, α3=3.3))

    @test s12_e + s3_e == s123_e
    @test s12_e + s3 == s123_e
    @test s12 + s3_e == s123_e
end

@testset "Multiplication: c * FSum, FSum * c " begin
    @test d == (α1=0.1,) * g1 + (α2=0.2,) * g2 + (α3=0.3,) * g3
    x2 = (Ext(α1=0.1,) * g1 + (α2=0.2,) * g2 + (α3=0.3,) * g3) |> d->fixpar(d,:α1)
    @test typeof(d2) == typeof(x2) && d2.fs == x2.fs && d2.αs == x2.αs
    @test nd == (α1=0.1,) * ng1 + (α2=0.2,) * ng2 + (α3=0.3,) * ng3    
end