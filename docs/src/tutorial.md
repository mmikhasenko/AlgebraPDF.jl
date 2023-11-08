```julia
using Markdown
using InteractiveUtils

begin
	using Parameters
	using AlgebraPDF
	using AlgebraPDF.QuadGK
    #    
    using Optim
	using Plots
	# 
	using Random
	Random.seed!(100)
end

md"""
# AlgebraPDF tutorial
"""

theme(:wong, frame=:box, lab="",
	xlim=(:auto,:auto), ylim=(0,:auto), grid=false)

md"""
## Gaussian
"""

gaussian = FGauss((μ=1.1, σ=0.9))

plot(gaussian, -4, 7, fill=0, α=0.2)

nGaussian = Normalized(gaussian, (-4, 7))

plot(nGaussian, fill=0, α=0.2)

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

model(0.0)

md"""
when array is passed, the function is broadcasted
"""

model(-1.8:0.9:1)

md"""
One can call it parameters
"""

model(0.0; p=(; α = -0.2, μ = 1.1, σ = 0.9, N1 = 0.85, N2 = 0.15))

md"""
The parameters can be adjusted
"""

model(0.0; p=(; α = -0.5, μ = 1.1, σ = 0.9, N1 = 0.85, N2 = 0.15))

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
