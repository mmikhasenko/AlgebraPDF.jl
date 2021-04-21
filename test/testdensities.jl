

#              _|_|_|                                
#    _|_|_|  _|          _|_|_|  _|    _|    _|_|_|  
#  _|    _|  _|  _|_|  _|    _|  _|    _|  _|_|      
#  _|    _|  _|    _|  _|    _|  _|    _|      _|_|  
#    _|_|_|    _|_|_|    _|_|_|    _|_|_|  _|_|_|    


@testset "Precodded pdfs" begin
    pdf1 = aGauss((mΩb = 6030, σ=17.0), (5600, 6400))
    @test freepars(pdf1) === (mΩb = 6030, σ=17.0)
    pdf2 = aBreitWigner((mΩb = 6030, Γ=17.0), (5600, 6400))
    @test freepars(pdf2) === (mΩb = 6030, Γ=17.0)
    pdf3 = aExp((τ = -1.1,), (-2, 2))
    @test freepars(pdf3) === (τ = -1.1,)

    pdfPowExp = aPowExp((n=3.3, τ=-4.0), (0, 2.2))
    @test pdfPowExp(1.1) != 0.0
    @test freepars(pdfPowExp) === (n=3.3, τ=-4.0)
    # 
    pars3 = (c0=1.1, c1=2.2, c2=3.3, c4=4.4)
    pdfPpol3 = aPol(pars3, (-2,12))
    @test func(pdfPpol3, 1.1; p=pars3) == sum(1.1^(i-1)*c for (i,c) in enumerate(pars3))
    @test freepars(pdfPpol3) === pars3
    # 
    pdf4 = aDoubleGaussFixedRatio((m = 0.77, Γ=0.15), (0, 1.0); fixpars=(r=0.8,n=3))
    @test pdf4(0.77) != 0.0
    pdf5 = aBreitWignerConvGauss((m = 0.77, Γ=0.15), (0, 1.0); fixpars=(σ=0.03,))
    @test pdf5(0.77) != 0.0
    # 
    xv = range(-π, 2π, length=40)
    yv = map(x->x*cos(3x) - 3*sin(x), xv)
    pdf6 = aTabulated(xv,yv,(-π,π))
    @test length(freepars(pdf6)) == 0
    @test pdf6(1.1) != 0.0
    @test pdf6(3π) == 0.0
    @test pdf6(-π) != 0.0
    @test pdf6(2π) != 0.0
end

# here is MWE of the normalization problem
# let 
#     d = aGauss((a=0.01,b = 1.1), (-3,3))
#     d(1.1, [36.1,0.1])
# end

