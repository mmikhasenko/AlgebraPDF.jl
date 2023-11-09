using AlgebraPDF, AlgebraPDF.Parameters

using Plots
theme(:wong, frame=:box, xlab="x", lab="", minorticks=true,
	guidefontvalign=:top, guidefonthalign=:right,
	xlim=(:auto,:auto), ylim=(0,:auto), grid=false)

myf(x; p=(a=1.1,b=2.2)) = x*p.a + p.b/x

gaussian = FGauss((μ=1.1, σ=0.9))

gaussian(1.1)

gaussian(-1.8:0.9:1)

pars(gaussian)

gaussian(0.0; p=(; μ = 1.1, σ = 0.9))

gaussian(0.0; p=(; μ = 1.1, σ = 0.9))

plot(gaussian, -4, 7, fill=0, α=0.8)
savefig("gaussian.pdf")

nGaussian = Normalized(gaussian, (-4, 7))

@assert nGaussian == gaussian |> Normalized((-4, 7))

nGaussian(1.1)

nGaussian(-1.8:0.9:1)

nGaussian(0.0; p=(; μ = 1.1, σ = 0.9))

plot(nGaussian, fill=0, α=0.7)
savefig("nGaussian.pdf")
