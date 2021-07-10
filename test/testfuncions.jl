using AlgebraPDF
using Test


@testset "FunctionWithParameters{NamedTuples}" begin
    d0 = FunctionWithParameters(
        (x;p)->p.a+cos(x)*p.b;
        p=(a=2,b=1))
    # 
    @test freepars(d0) == (a=2,b=1)
    @test fixedpars(d0) == ∅
    @test pars(d0) == (a=2,b=1)

    @test nfreepars(d0) == 2

    @test func(d0,0) == 3
    @test func(d0,repeat([0],10)) == repeat([3],10)
    @test func(d0, 1:10) == func(d0, collect(1:10))

    d1 = updatepar(d0, :a, 4)
    @test freepars(d1) == (a=4,b=1)
    @test releasepar(d0, :b) == d0
    @test_throws ArgumentError fixpar(d0, :b)
    # 
    @test v2p([1,2],d0) == (a=1,b=2)
    @test p2v(d0) == [2,1]
    @test p2v((b=6,a=5,c=33), d0) == [5,6]    
end

@testset "FunctionWithParameters{FlaggedNamedTuple}" begin
    g = FunctionWithParameters(
        (x;p)->exp((x-p.μ)^3/(2*p.σ^2));
        p=Ext(μ=1.2, σ=0.1))
    # 
    g′ = fixpar(g, :μ)
    @test fixedpars(g′).μ == 1.2
    g′ = fixpar(g, :μ, 1.3)
    @test fixedpars(g′).μ == 1.3
    # 
    @test freepars(g′) == (σ=0.1,)
    @test fixedpars(g′) == (μ=1.3,)
    #
    @test g′(1.3; p=(σ=0.2,)) == 1.0 # now works since μ is fixed
    @test g′(1.3; p=(σ=0.2, μ=1.6)) == 1.0 # will use the fixed value of μ=1.2
    g′′ = updatepar(g′, :μ, 1.2) # however, the update fill do
    @test g′′(1.2) == 1.0
    # 
    g′′′ = releasepar(g′′, :μ)
    @test g′′′ == g # true
end

@testset "FunctionWithParameters{FlaggedNamedTuple}" begin
    m0 = FunctionWithParameters(
        (x;p)->p.a+cos(x)*p.b;
        p=Ext(a=2,b=1))
    #
    m1 = updatepar(m0, :a, 4)
    @test freepars(m1) == (a=4,b=1)

    m2 = fixpar(m1, :a)
    @test freepars(m2) == (b=1,)
    @test fixedpars(m2) == (a=4,)
    # 
    m3 = releasepar(m0, :a)
    @test m3 == m0
    # 
    fixpars(m1, (:a,:b)) == fixpars(m1, [:a,:b])
    # 
    @test m1 == releasepars(m2, (:a,))
    @test m1 == releasepars(m2, [:a])
    # 
    m4 = updatepars(m0, (a=4, b=5))
    @test freepars(m4) == (a=4,b=5)
end

@testset "func on scalars" begin
    scalar_function = 1.1
    @test func(scalar_function, 3.3) == 1.1
    @test pars(scalar_function, rand([true,false])) == ∅
    # 
    Function_function = sin
    @test func(Function_function, 0.0) ≈ 0.0
    @test pars(Function_function, rand([true,false])) == ∅
end


@testset "Abs2Func" begin
    m0 = FunctionWithParameters(
        (x;p)->p.a+cos(x)*p.b;
        p=Ext(a=2,b=1))
    # 
    am0 = abs2(m0)
    @test func(am0,π/2) == 4 
    am1 = updatepar(am0, :a, 5)
    @test func(am1,π/2) == 25
    @test freepars(am0) == freepars(m0)
end

@testset "LogFunc" begin
    m0 = FunctionWithParameters(
        (x;p)->exp(p.b*x);
        p=Ext(b=2,))
    # 
    am0 = log(m0)
    @test am0(1) == 2 
    am1 = updatepar(am0, :b, 4)
    @test am1(π/2) ≈ 2π
    @test freepars(am0) == freepars(m0)
end

@testset "NegativeLogLikelihood" begin
    d = Normalized(FGauss((μ=1.1, σ=0.3)), (-3, 5))
    data = randn(1000)
    nll = NegativeLogLikelihood(d, data)
    # 
    @test nll == minussum(log(d), data)
    @test nll(0.0) != 0.0
    @test nll(0.0) == nll(10.0)
    @test nll(0.0; p=(μ=1, σ=2)) == nll(0.0, [1,2])
end

@testset "SumFunc and ProdFunc" begin
    a1 = FunctionWithParameters(
        (x;p)->p.a+cos(x)*p.b; p=Ext(a=2,b=1))
    a2 = FunctionWithParameters(
        (x;p)->p.a+sin(x)*p.b; p=Ext(a=2,b=1))
    #
    sum1 = a1+a2
    sum2 = +(a1,a2,(c1=1.0,c2=1.0))
    sum3 = +(a1,a2,Ext(d1=1.0, d2=1.0))

    @test func(sum1,2) == func(sum2,2)
    @test func(sum1,2) == func(sum3,2)

    @test freepars(sum1) == (a=2, b=1, α1=1.0, α2=1.0)
    @test freepars(sum2) == (a=2, b=1, c1=1.0, c2=1.0)
    @test freepars(sum3) == (a=2, b=1, d1=1.0, d2=1.0)
    # 
    @test freepars(fixpar(sum3,:a)) == (b=1, d1=1.0, d2=1.0)
    # 
    sum4 = fixpar(sum3, :a, 3)
    @test fixedpars(sum4) == (a=3,)
    @test fixedpars(sum4[1]) == (a=3,)
    @test fixedpars(sum4[2]) == (a=3,)
    #
    sum5 = fixpar(sum3,:d2, 2)
    @test func(sum5,3.3) == func(a1,3.3) + 2* func(a2,3.3)
    # 
    sum6 = a1-a2
    @test sum6(1.1) == a1(1.1)-a2(1.1)
    # 
    prd = a1*a2
    @test prd(1.1) == a1(1.1)*a2(1.1)
end


@testset "no-parameters f and no-parameters normalized f" begin
    d0 = FunctionWithParameters((e;p)->e^2+p.a; p=(a=1.0,))
    f0 = noparsf(d0)
    xv = -1:10
    @test func(d0,xv;p=freepars(d0)) ≈ f0(xv)
end


# implementation with NAMES of parameters build into the funciton call
struct BW1{P} <: AbstractFunctionWithParameters
    p::P
end
import AlgebraPDF:func
func(bw::BW1, x::NumberOrTuple; p=pars(bw)) = p.m*p.Γ/(p.m^2-x^2-1im*p.m*p.Γ)


@testset "User-def type with fixed name" begin
    bw = BW1((m=3.1,Γ=0.1))
    @test pars(bw).m == 3.1
    @test pars(bw).Γ == 0.1
    @test bw(1.1) != 0.0
end


# implementation with ORDER of parameters build into the funciton call
struct BW2{P} <: AbstractFunctionWithParameters
    p::P
end
import AlgebraPDF:func
function func(bw::BW2, x::NumberOrTuple; p=pars(bw))
    m,Γ = (getproperty(p,s) for s in keys(bw.p))
    m*Γ/(m^2-x^2-1im*m*Γ)
end


@testset "User-def type with fixed order" begin
    bw = BW2((m1=3.1,Γ1=0.1))
    # 
    @test pars(bw).m1 == 3.1
    @test pars(bw).Γ1 == 0.1
    @test real(bw(3.1)) ≈ 0.0
    @test real(bw(3.1; p=(m1=1.1,Γ1=3.3))) != 0
    # 
    bw = BW2((m2=3.1,Γ2=0.1))
    @test pars(bw).m2 == 3.1
    @test pars(bw).Γ2 == 0.1
end

@makefuntype SuperF(x;p) = x^2+p.a*x^3
g = SuperF(p=(a=0.5,))

@testset "@maketype macro" begin
    @test keys(pars(g)) == (:a,)
    @test g(1) == 1.5
end
