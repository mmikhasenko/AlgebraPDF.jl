using AlgebraPDF, AlgebraPDF.Parameters
using LinearAlgebra, Optim

using Random
Random.seed!(100)

using Plots
theme(:wong, frame=:box, xlab="x", lab="", minorticks=true,
	guidefontvalign=:top, guidefonthalign=:right,
	xlim=(:auto,:auto), ylim=(0,:auto), grid=false)

myf(x; p=(a=1.1,b=2.2)) = x*p.a + p.b/x ;

gaussian = FGauss((μ=1.1, σ=0.9))

gaussian(1.1)

gaussian(-1.8:0.9:1)

pars(gaussian)

gaussian(0.0; p=(; μ = 1.1, σ = 0.9))

gaussian(0.0; p=(; μ = 0.0, σ = 1.9))

plot(gaussian, -4, 7, fill=0, α=0.8)
savefig("gaussian.pdf")

nGaussian = Normalized(gaussian, (-4, 7))

@assert nGaussian == gaussian |> Normalized((-4, 7))

nGaussian(1.1)

nGaussian(-1.8:0.9:1)

nGaussian(0.0; p=(; μ = 1.1, σ = 0.9))

plot(nGaussian, fill=0, α=0.7)
savefig("nGaussian.pdf")

exponential = FExp((; α=-0.2))

nExponent = exponential |> Normalized((-4, 7))

plot(nExponent)
savefig("nExponential.pdf")

model = FSum([nExponent, nGaussian], (N1=0.85, N2=0.15))

@assert model == nExponent * (N1=0.85, ) + nGaussian * (N2=0.15,)

begin
    plot(model)
    plot!(model[1], ls=:dash, lab="background")
    plot!(model[2], fill=0, lab="signal")
end
savefig("model_components.pdf")

@time data = AlgebraPDF.rand(model, 10_000)

stephist(data, bins=100)
savefig("sampled_data_histogram.pdf")

nll = NegativeLogLikelihood(model, data)

nll(1.0), nll(0.0), nll(())

ext = Extended(nll)

starting_values = let
    Nd = length(data)
    default_values = pars(ext)
    @unpack N1, N2 = default_values
    Nsum = N1+N2
    merge(default_values, (N1 = N1/Nsum*Nd, N2 = N2/Nsum*Nd))
end

@time fit = let
    initial_invH = Diagonal([0.001,0.01,0.01,100,100]) .+ eps()

    optimize(x->ext(1.1, x), starting_values |> collect,
	    BFGS(; initial_invH = x -> initial_invH,))
end

best_model = updatepars(model, NamedTuple{keys(pars(model))}(fit.minimizer))

let
    bins = range(lims(model)..., 100)
    Nd = length(data)

    stephist(data; bins)
    plot!(model, scaletobinneddata(Nd, bins), lab="original")
    plot!(best_model, scaletobinneddata(bins), lab="fit")

    plot!(best_model[2], scaletobinneddata(bins), fill=0, lab="signal")
    plot!(best_model[1], scaletobinneddata(bins), ls=:dash, lab="background")
end
savefig("model_fitting.pdf")
