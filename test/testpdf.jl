using AlgebraPDF
using Test

@testset "Construction of pdf with parametric function" begin
    d0 = Normalized(FunctionWithParameters((x;p)->x^2+p.a; p=(a=2.0,)), (-2, 4))
    func(d0, 3) == 11
    pars(d0) == (a=2.0,)
    d1 = updatepars(d0, (a=3.0,))
    func(d1, 3) == 12
end

@testset "short cut: Normalized(function), fixedshapepdf" begin
    d = Normalized((e;p)->e^2+p.a; lims=(-1,2), p=(a=1.0,))
    @test typeof(pars(d)) <: NamedTuple
    @test length(freepars(d)) == 1
    @test length(fixedpars(d)) == 0
    # 
    d1 = fixedshapepdf(x->exp(-(4x)^2), (-1, 2))
    @test length(freepars(d1)) == 0
end

@testset "fix parameters example" begin
    d0 = Normalized((e;p)->e^2+p.a; p=Ext(a=1.0), lims=(-1,2))
    @test d0(1.0) ≈ 2.0/(9/3+3*pars(d0).a)
    #
    d1 = fixpar(d0, :a, 2.9)
    @test d1(1.0) ≈ 3.9/(9/3+3*pars(d1).a)
    # 
    d1′ = fixpars(d0, (:a,))
    @test freepars(d1′) == freepars(fixpars(d0, [:a]))
    @test freepars(d1′) == freepars(fixpars(d0, (a=1.0,)))
end

@testset "noparsnormf " begin
    d = Normalized((e;p)->e^2+p.a; p=Ext(a=1.0), lims=(-1,2))
    f = noparsnormf(d)
    d(1.1) == f(1.1)
end

@testset "integral " begin
    d1 = Normalized((e;p)->e^2+p.a; p=(a=1.0,), lims=(-1,2))
    @test integral(d1, (0,2)) < 1
    d2 = Normalized((e;p)->-e^2+p.b; p=(b=4.0,), lims=(-1,2))
    @test integral(d2, lims(d2)) ≈ 1
    s = d1+d2
    @test integral(s) == 2
    @test 1.8 < integral(s, (-0.9,1.9)) < 2
    s2 = updatepars(s,(α1=2,α2=4))
    # 
    @test integral(s2[1]) == 2
    @test integral(s2[2]) == 4
    # 
    @test 1.7 < integral(s2[1], (-0.9,1.9)) < 2
end

# implementation with NAMES of parameters build into the funciton call
struct nBW1{P,L} <: AbstractPDF{1}
    p::P
    lims::L
end
import AlgebraPDF:func
func(bw::nBW1, x::NumberOrTuple; p=pars(bw)) = abs2(p.m*p.Γ/(p.m^2-x^2-1im*p.m*p.Γ))

@testset "User-def type with fixed name" begin
    bw = nBW1((m=3.1,Γ=0.1), (0,4))
    @test pars(bw).m == 3.1
    @test pars(bw).Γ == 0.1
    @test bw(3.1) != 1.0 # because of normalization
    @test lims(bw) == (0,4)
end


# implementation with ORDER of parameters build into the funciton call
struct nBW2{P,L} <: AbstractPDF{1}
    p::P
    lims::L
end
import AlgebraPDF:func
function func(bw::nBW2, x::NumberOrTuple; p=pars(bw))
    m,Γ = (getproperty(p,s) for s in keys(bw.p))
    m*Γ/(m^2-x^2-1im*m*Γ)
end


@testset "User-def type with fixed order" begin
    bw_i = nBW2((m_i=3.1,Γ_i=0.1), (0.0,3.0))
    @test pars(bw_i).m_i == 3.1
    @test pars(bw_i).Γ_i == 0.1
    @test lims(bw_i) == (0.0,3.0)
    bw_j = nBW2((m_j=4.1,Γ_j=0.2), (1,4))
    @test pars(bw_j).m_j == 4.1
    @test pars(bw_j).Γ_j == 0.2
    @test lims(bw_j) == (1,4)
end

# same as the first option but with a macro
@makepdftype nBW3(y;p) = abs2(p.m*p.Γ/(p.m^2-y^2-1im*p.m*p.Γ))

@testset "User-def type with fixed name" begin
    bw = nBW3((m=3.1,Γ=0.1), (0,4))
    @test pars(bw).m == 3.1
    @test pars(bw).Γ == 0.1
    @test bw(3.1) != 1.0 # because of normalization
    @test lims(bw) == (0,4)
end
