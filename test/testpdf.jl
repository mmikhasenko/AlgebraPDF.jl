
@testset "Parameters of pdf" begin
    d = pdf((e;p)->e^2+p.a; lims=(-1,2), p=(a=1.0,))
    @test typeof(pars(d)) <: NamedTuple
    @test length(freepars(d)) == 1
    @test length(fixedpars(d)) == 0
end

@testset "fix parameters example" begin
    d0 = pdf((e;p)->e^2+p.a; p=TwoNamedTuples(a=1.0), lims=(-1,2))
    @test d0(1.0) ≈ 2.0/(9/3+3*pars(d0).a)
    #
    d1 = fixpar(d0, :a, 2.9)
    @test d1(1.0) ≈ 3.9/(9/3+3*pars(d1).a)
    # 
    d1′ = fixpars(d0, (:a,))
    @test freepars(d1′) == freepars(fixpars(d0, [:a]))
    @test freepars(d1′) == freepars(fixpars(d0, (a=1.0,)))
end

@testset "parameters to values" begin
    d = pdf((e;p)->e^2+p.a; p=(a=1.0,), lims=(-1,2))
    @test p2v(d) == [1.0]
    @test p2v((a=3.0,), d) == [3.0]
end

@testset "fixed-shape pdf" begin
    d1 = fixedshapepdf(x->exp(-(4x).^2), (-1, 2))
    @test length(freepars(d1)) == 0
end

@testset "no-parameters f and no-parameters normalized f" begin
    d0 = pdf((e;p)->e^2+p.a; p=(a=1.0,), lims=(-1,2))

    f = noparsf(d0)
    xr = lims(d0)[1]+rand()*(lims(d0)[2]-lims(d0)[1])
    @test func(d0,xr;p=freepars(d0)) ≈ f(xr)
    #
    ananorm = ((8+1)/3+1.0*3)
    @test d0(xr) ≈ f(xr)/ananorm
end


# @testset "example of usage" begin
#     # tests
#     snl = pdf((x;p) -> exp(-(x-p.μ)^2/(2*p.σ^2)); p = (μ=1.4, σ=0.15), lims=(0, 3))
#     bkg = pdf((x;p) -> sqrt(x)*exp(-p.α*x); p = (α=1.3,), lims=(0, 3))
#     bkg *= (fb=2.5,)
#     pdf_sum = snl + bkg
    
#     # generating
#     data = generate(10000, pdf_sum; p=freepars(pdf_sum));
#     @test length(data) == 10000
    
#     # fitting
#     fr = fit_llh(data, pdf_sum)
#     pfr = v2p(minimizer(fr), pdf_sum)
#     @test length(pfr) == 4
#     # 
#     @test minimum(fr) ≤ llh(data, pdf_sum)
#     # 
#     invH_bfgs = invH(fr)
#     @test size(invH_bfgs) == (4,4)
#     @test covmat(fr) == invH_bfgs
#     @test length(errors(fr)) == 4
#     # 
#     invH_fd = invexacthessian(fr)
#     rel_error_max = max((abs.(diag(invH_bfgs) .- diag(invH_fd)) ./ abs.(diag(invH_fd)))...)
#     @test rel_error_max < 10e-2 # 10%
#     # 
#     mfr = measurements(fr)
#     @test typeof([measurements(fr)...]) <: Vector{Measurement{T}} where T
#     mfr_exact = measurements(fr, true)
#     @test Measurements.uncertainty.(mfr) == errors(fr)
#     @test Measurements.uncertainty.(mfr_exact) == sqrt.(diag(invH_fd))
# end

@testset "pdf with NamedTuple" begin
    # here I create the structure directly,
    #    p=(...) would call TwoNamedTuples(p)
    d = FunctionWithParameters((x;p)->x*5+p.a, (a=1.1,))
    # 
    @test func(d,1.1) != 0.0
    #
    # parameter manipulation are not supposed to work with NamedTuple
    @test_throws ArgumentError fixpar(d, :a, 1.2)
    @test_throws ArgumentError fixpars(d, (a=1.2,))
end


#  _|_|_|  _|_|      _|_|_|    _|_|_|  _|  _|_|    _|_|      _|_|_|  
#  _|    _|    _|  _|    _|  _|        _|_|      _|    _|  _|_|      
#  _|    _|    _|  _|    _|  _|        _|        _|    _|      _|_|  
#  _|    _|    _|    _|_|_|    _|_|_|  _|          _|_|    _|_|_|    


@macroexpand( @typepdf BW(x; p) = p.m*p.Γ/(p.m^2-x^2-1im*p.m*p.Γ) )


@macroexpand( @typepdf BW(x; p) = p.m*p.Γ/(p.m^2-x^2-1im*p.m*p.Γ) )

@typepdf BW(x; p) = p.m*p.Γ/(p.m^2-x^2-1im*p.m*p.Γ) |> abs2
bw = BW((m=3.1,Γ=0.1), (1,5))

@testset "User-def type" begin
    @test pars(bw).m == 3.1
    @test pars(bw).Γ == 0.1
    @test lims(bw) == (1,5)
    @test normalizationintegral(bw) != 0.0
    @test bw(1.1) != 0.0
end
