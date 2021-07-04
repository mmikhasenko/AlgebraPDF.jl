using SpecialFunctions

#    _|_|_|    _|_|    _|_|_|    _|      _|  
#  _|        _|    _|  _|    _|  _|      _|  
#  _|        _|    _|  _|    _|    _|  _|    
#    _|_|_|    _|_|    _|    _|      _|      

@makepdftype myBW(x;p) = let xth = 2.961477
    AlgebraPDF.amplitudeBWsq(x+xth, p.Δm+xth, p.Γ)
end
mybw = myBW(Ext(Δm=0.1, Γ=3e-3), (0,0.22))

# convolute with average resolution
const σ_mean = 5e-3
mybw_conv = convGauss(mybw, σ_mean)

@testset "convGauss: update pars works" begin
    @test pars(updatepars(mybw_conv, (Δm=0.15,))).Δm == 0.15
end

function updategetmaxdensity(d, Δm)
    xv = range(0,0.22,length=100)
    maximum(updatepars(d,(Δm=Δm,))(xv))
end

# resolution is linear function of energy
resolution_model = FunctionWithParameters((x;p)->σ_mean+(x-0.1)*3e-2; p=∅)
mybw_conv_dep = convGauss(mybw, resolution_model)

@testset "convGauss energy-dependent gets wider" begin
    # if the peak located below 0.1 the resolution the effective resolution  < σ_mean => the peak max is higher
    @test updategetmaxdensity(mybw_conv_dep, 0.05) > updategetmaxdensity(mybw_conv, 0.05)
    # if the peak located above 0.1 the resolution the effective resolution  > σ_mean => the peak max is smaller
    @test updategetmaxdensity(mybw_conv_dep, 0.15) < updategetmaxdensity(mybw_conv, 0.15)
end
