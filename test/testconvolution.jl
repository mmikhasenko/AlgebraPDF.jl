#    _|_|_|    _|_|    _|_|_|    _|      _|  
#  _|        _|    _|  _|    _|  _|      _|  
#  _|        _|    _|  _|    _|    _|  _|    
#    _|_|_|    _|_|    _|    _|      _|      


function df(f1,f2,lims; Ns = 10)
    xv = range(lims..., length=Ns)
    dyv = f1.(xv) .- f2.(xv)
    return sqrt(sum(abs2, dyv) / Ns)
end


@testset "convolution with gauss" begin
    @test AlgebraPDF.standardgauss(0,1) ≈ 1/sqrt(2π)
    @test AlgebraPDF.standardgauss(2,2) ≈ exp(-1/2)/sqrt(2π*4)
    # tests
    σ0 = 0.3
    aconv(e) = (erf(e/sqrt(2*σ0^2))+1)/2
    #
    nconv1(e) = conv_with_gauss(e, x->x>0, σ0)
    @test df(nconv1, aconv, (-1, 1); Ns=1000) < 0.01
    #
    step = pdf((e;p)->e>0; p=∅, lims=(-1,1))
    func(step, 1.1; p=∅)
    smeared_step = conv_with_gauss(step, σ0)
    nconv2(e) = func(smeared_step, e; p=∅)
    # smeared_step(1.1)
    func(smeared_step, 1.1; p=∅)
    @test df(nconv2, aconv, (-1, 1); Ns=1000) < 0.01
    # 
    smeared_step_sampling = conv_with_gauss_sampling(step, σ0; Ns=50)
    nconv3(e) = func(smeared_step_sampling, e; p = freepars(smeared_step_sampling))
    @test df(nconv3, aconv, (-1, 1); Ns=1000) < 0.01
    # 
    smeared_step_sampling = conv_with_gauss_sampling(step, σ0; Ns=10)
    nconv4(e) = func(smeared_step_sampling, e; p = freepars(smeared_step_sampling))
    @test 0.01 < df(nconv4, aconv, (-1, 1); Ns=1000) < 0.04
end


@typepdf myBW(x;p) = AlgebraPDF.amplitudeBWsq(x+2.961477, p.Δm+2.961477, p.Γ)
mybw = myBW(TwoNamedTuples(Δm=0.1, Γ=3e-3), (0,0.22))
mybw_conv = convGauss(mybw, 5e-3)

@testset "convGauss: update pars works" begin
    @test pars(updatepars(mybw_conv, (Δm=0.15,))).Δm == 0.15
end


function maxdensity(d, Δm)
    xv = range(0,0.22,length=100)
    maximum(updatepars(d,(Δm=Δm,))(xv))
end

resolution_model = FunctionWithParameters((x;p)->5e-3+(x-0.1)*3e-2; p=∅)
mybw_conv_dep = convGauss(mybw, resolution_model)
@testset "convGauss energy dep. gets wider" begin
    @test maxdensity(mybw_conv_dep, 0.05) > maxdensity(mybw_conv, 0.05)
    @test maxdensity(mybw_conv_dep, 0.15) < maxdensity(mybw_conv, 0.15)
end