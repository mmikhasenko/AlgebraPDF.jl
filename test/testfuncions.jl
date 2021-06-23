using AlgebraPDF

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

@testset "FunctionWithParameters{TwoNamedTuples}" begin
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

@testset "SumFunc" begin
    a1 = FunctionWithParameters(
        (x;p)->p.a+cos(x)*p.b; p=Ext(a=2,b=1))
    a2 = FunctionWithParameters(
        (x;p)->p.a+sin(x)*p.b; p=Ext(a=2,b=1))
    #
    sum1 = a1+a2
    sum2 = +(a1,a2,(c=1.0,))
    sum3 = +(a1,a2,Ext(d=1.0))

    @test func(sum1,2) == func(sum2,2)
    @test func(sum1,2) == func(sum3,2)

    @test freepars(sum1) == (a=2, b=1, α2=1.0)
    @test freepars(sum2) == (a=2, b=1, c=1.0)
    @test freepars(sum3) == (a=2, b=1, d=1.0)
    # 
    @test freepars(fixpar(sum3,:a)) == (b=1, d=1.0)
    # 
    sum4 = fixpar(sum3,:a, 3)
    @test fixedpars(sum4) == (a=3,)
    @test fixedpars(sum4.f1) == (a=3,)
    @test fixedpars(sum4.f2) == (a=3,)
    # 
    sum5 = fixpar(sum3,:d, 2)
    @test func(sum5,3.3) == func(a1,3.3) + 2* func(a2,3.3)
end
