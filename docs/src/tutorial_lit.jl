#md # # Tutorial
#nb # # AlgebraPDF.jl Tutorial
#jl # # AlgebraPDF.jl Tutorial

# This tutorial guides the user through the basic example of
# creating a density function which is a sum of a gaussian signal peak,
# and exponential background.

using AlgebraPDF, AlgebraPDF.Parameters

using Plots
theme(:wong, frame=:box, xlab="x", lab="", minorticks=true, 
	guidefontvalign=:top, guidefonthalign=:right,
	xlim=(:auto,:auto), ylim=(0,:auto), grid=false)

# ## Function With Parameters

# The package provides a standard wrapper of a function with parameters `f(x;p)`,
# where `x` is function variable, `p` is a structure that holds parameters
# with their names. The simplest and the most common case is where `x` would be a number, 
# and `p` is a named tuple. For example,

myf(x; p=(a=1.1,b=2.2)) = x*p.a + p.b/x

# The modele introduces a type `FunctionWithParameters`,
# which intends to behave like the `myf` from the user prospective.

# ## A Gaussian function

# The gaussian function is constructed calling a specific type `FGauss`,
# and giving a tuple of parameters with their default values 

gaussian = FGauss((μ=1.1, σ=0.9))

# It is on of the predefined examples of functions with parameters, `FGauss <: AbstractFunctionWithParameters`. 
# The object is callable like a regular function

gaussian(1.1)

# when array is passed, the function is broadcasted

gaussian(-1.8:0.9:1)

# Default values of the parameters can be accessed with `pars`, and `freepars` method.

pars(gaussian)

# One can provide the key argument `p` with the named tuple of parameters.
# These tuple is always used instread of the default values.

gaussian(0.0; p=(; μ = 1.1, σ = 0.9))

# The parameters can be adjusted

gaussian(0.0; p=(; μ = 1.1, σ = 0.9))

# Similar to the regular fuction, the object can be plotted.

plot(gaussian, -4, 7, fill=0, α=0.8)
#jl savefig("gaussian.pdf")
#md savefig("gaussian.pdf")
#md savefig("gaussian.svg"); nothing # hide
#md # [![gaussian](gaussian.svg)](gaussian.pdf)


# ## Normalization

# To turn an arbitrary function to the probability density,
# one need to introduce normalization.
# This is done by attaching the range (support) to the function.

nGaussian = Normalized(gaussian, (-4, 7))

# The same is archived with the pipeline.

@assert nGaussian == gaussian |> Normalized((-4, 7))

# The normalized object has the call method regular function.
# The normalization is computed on fly.
# It is constly for a single-value call.

nGaussian(1.1)

#  When calling on iterable collection, the normalization is computed once.

nGaussian(-1.8:0.9:1)

# As before, the parameters can be updates by passing a named tuple

nGaussian(0.0; p=(; μ = 1.1, σ = 0.9))

# For plotting of the normalized function, one does not need to specify the range.

plot(nGaussian, fill=0, α=0.7)
#jl savefig("nGaussian.pdf")
#md savefig("nGaussian.pdf")
#md savefig("nGaussian.svg"); nothing # hide
#md # [![nGaussian](nGaussian.svg)](nGaussian.pdf)

