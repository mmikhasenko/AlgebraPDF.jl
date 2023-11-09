using Optim
using Random
Random.seed!(100)


md"""
## Exponent
"""

background = FExp((; α=-0.2))

nExponent = background |> Normalized((-4, 7))

plot(nExponent)

md"""
## Sum
"""

model = FSum([nExponent, nGaussian], (N1=0.85, N2=0.15))

begin
	plot(model)
	plot!(model[1], ls=:dash, lab="background")
	plot!(model[2], fill=0, lab="signal")
end

md"""
## Parameters
"""

model # let's explore

md"""
Here are the default values of the parameters
"""

pars(model)


md"""
The model is immutable, but we can create a new model with different parameters
"""

begin
	plot()
	plot!(model, lab="default")
	plot!(updatepar(model, :μ, 1.4), lab="updated")
	plot!(updatepars(model, (μ=1.4, N1=0.5)), lab="much updated")
end

md"""
## Sampling
"""

md"""
Sampling is implemented by numerical Inversion Method
"""

@time data = AlgebraPDF.rand(model, 10_000)

stephist(data, bins=100)

md"""
The grid size can be ajusted with the `generate` method
"""

let Nsamples = 100_000
	df = DataFrame(Nbins = [10, 100, 1000])
	df.sampling_time = map(df.Nbins) do Nbins
		@elapsed generate(model, Nsamples; Nbins)
	end
	df
end

ρ = AlgebraPDF.getbinned1dDensity(x->model(x), lims(model), 10)

plot(ρ.grid, ρ.cumarr, m=(4, :o))

rand(ρ)

stephist([rand(ρ) for _ in 1:100_000], bins=100)

md"""
## Likelihood
"""

nll = NegativeLogLikelihood(model, data)

ext = Extended(nll)

md"""
The size of the data 10_000. Also, it is normalized
"""

@time nll(())

md"""
The negative log likelihood does not depend on x, only depend on parameters
"""

nll(1.0), nll(0.0), nll(())

starting_values = let
	Nd = length(data)
	default_values = pars(ext)
	@unpack N1, N2 = default_values
	Nsum = N1+N2
	# 
	default_values + (N1 = N1/Nsum*Nd, N2 = N2/Nsum*Nd)
end

ext(1.1, collect(starting_values)+[0,0,0,0,0])

@time fit = let
	initial_invH = Diagonal([0.001,0.01,0.01,100,100]) .+ eps()
	
	optimize(x->ext(1.1, x), starting_values |> collect,
	BFGS(; initial_invH = x -> initial_invH,))
end

best_model = updatepars(model, NamedTuple{keys(pars(model))}(fit.minimizer))

let
	bins = range(lims(model)..., 100)
	Nd = length(data)
	# 
	stephist(data; bins)
	plot!(model, scaletobinneddata(Nd, bins), lab="original")
	plot!(best_model, scaletobinneddata(bins), lab="fit")
	# 
	plot!(best_model[2], scaletobinneddata(bins), fill=0, lab="signal")
	plot!(best_model[1], scaletobinneddata(bins), ls=:dash, lab="background")
end
```
