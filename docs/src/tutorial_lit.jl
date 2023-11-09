
#md This tutorial guides the user through the basic example of
#md creating a density function which is a sum of a gaussian signal peak,
#md and exponential background.

using AlgebraPDF
using AlgebraPDF.QuadGK
#    
using Parameters
using Plots


theme(:wong, frame=:box, lab="",
	xlim=(:auto,:auto), ylim=(0,:auto), grid=false)


#md The gaussian function is constructed calling a specific type `FGauss`,
#md and giving a tuple of parameters with their default values 

gaussian = FGauss((μ=1.1, σ=0.9))

#md The object is callable like a regular function

gaussian(1.1)

#md when array is passed, the function is broadcasted

gaussian(-1.8:0.9:1)

#md Default values of the parameters can be accessed with `pars`, and `freepars` method.

pars(gaussian)
freepars(gaussian)

#md One can provide the key argument `p` with the named tuple of parameters.
#md These tuple is always used instread of the default values.

gaussian(0.0; p=(; μ = 1.1, σ = 0.9))

#md The parameters can be adjusted

gaussian(0.0; p=(; μ = 1.1, σ = 0.9))

#md Similar to the regular fuction, the object can be plotted.

plot(gaussian, -4, 7, fill=0, α=0.2)
#jl savefig("gaussian.pdf")
#md savefig("gaussian.pdf")
#md savefig("gaussian.svg"); nothing # hide
#md # [![Data](gaussian.svg)](gaussian.pdf)


#md ## Normalization

#md To turn an arbitrary function to the probability density,
#md one need to add normalization.
#md This is done by attaching the range (support) to the function.

nGaussian = Normalized(gaussian, (-4, 7))

#md The same is archived with the pipeline.

@assert nGaussian == gaussian |> Normalized((-4, 7))

#md The normalized object has the call method regular function.
#md The normalization is computed on fly.
#md It is constly for a single-value call.

nGaussian(1.1)

#md When calling on iterable collection, the normalization is computed once.

nGaussian(-1.8:0.9:1)

# As before, the parameters can be updates by passing a named tuple

nGaussian(0.0; p=(; μ = 1.1, σ = 0.9))


#md ## Plotting

#md For plotting of the normalized function, one does not need to specify the range.

plot(nGaussian, fill=0, α=0.7)
#jl savefig("nGaussian.pdf")
#md savefig("nGaussian.pdf")
#md savefig("nGaussian.svg"); nothing # hide
#md # [![Data and True Statistical Model](nGaussian.svg)](nGaussian.pdf)

