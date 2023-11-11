#md # # Tutorial
#nb # # AlgebraPDF.jl Tutorial
#jl # # AlgebraPDF.jl Tutorial

# This tutorial guides the user through the basic example of
# creating a density function which is a sum of a gaussian signal peak,
# and exponential background, sampling from the distribution,
# and optimizating its parameters with an unbinned fit.

using AlgebraPDF, AlgebraPDF.Parameters
using LinearAlgebra, Optim

using Plots
theme(:wong, frame=:box, xlab="x", lab="", minorticks=true, 
	guidefontvalign=:top, guidefonthalign=:right,
	xlim=(:auto,:auto), ylim=(0,:auto), grid=false)

# ## Function With Parameters

# The package provides a standard wrapper of a function with parameters `f(x;p)`,
# where `x` is function variable, `p` is a structure that holds parameters
# with their names. The simplest and the most common case is where `x` would be a number, 
# and `p` is a named tuple. For example,

myf(x; p=(a=1.1,b=2.2)) = x*p.a + p.b/x ;

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

gaussian(0.0; p=(; μ = 0.0, σ = 1.9))

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


# This section introduces the exponential function in `AlgebraPDF.jl` and demonstrates how to normalize it.
# Understanding how to work with different types of functions and their normalization 
# is crucial in creating complex distributions.

# Exponential function is another lineshape defined in the package.
# We are using the expnential distribution with a slope `α` to define the background to our gaussian signal.

exponential = FExp((; α=-0.2))

# Next, we normalize this function over a specific range in the same way we did with the Gaussian function in the previous sections.
# Normalization is essential to ensure that the function behaves as a probability density function
# and we can countably measure the background fraction

nExponent = exponential |> Normalized((-4, 7))

# Let's plot the normalized exponential function to visualize it.

plot(nExponent)
#jl savefig("nExponential.pdf")
#md savefig("nExponential.pdf")
#md savefig("nExponential.svg"); nothing # hide
#md # [![nExponential](nExponential.svg)](nExponential.pdf)

# ## Summation of Functions

# Now, we can explore how to sum different types of functions or distributions.
# Summing functions is essential part of the package functionality that allows one to model more complex distributions.

# Here, we create a model that is a sum of the previously defined normalized exponential function and the normalized Gaussian function.

model = FSum([nExponent, nGaussian], (N1=0.85, N2=0.15))

# The sum of functions, is an object of type `FSum`, that hold a list of functions and their weights in static vectors of equal sizes.
# There is an alternative way to formulate an equivalent model

@assert model == nExponent * (N1=0.85, ) + nGaussian * (N2=0.15,)

# where the product of the function with a named tuple return `FSum` object. The summation between two `FSum` objects is defined.
# Complementary, individual components with their weight can be accessed by indexing the `FSum` object as in the following protting code.

begin
    plot(model)
    plot!(model[1], ls=:dash, lab="background")
    plot!(model[2], fill=0, lab="signal")
end
#jl savefig("model_components.pdf")
#md savefig("model_components.pdf")
#md savefig("model_components.svg"); nothing # hide
#md # [![Model and Components](model_components.svg)](model_components.pdf)




# ## Sampling from the model

# Sampling from a model is a key operation in statistical modeling.
# `AlgebraPDF.jl` implements sampling through the numerical Inversion Method.

# Here, we sample data points from our model.

@time data = AlgebraPDF.rand(model, 10_000)

# Plotting the sampled data as a histogram gives us a visual representation of the distribution.

stephist(data, bins=100)
#jl savefig("sampled_data_histogram.pdf")
#md savefig("sampled_data_histogram.pdf")
#md savefig("sampled_data_histogram.svg"); nothing # hide
#md # [![Sampled Data Histogram](sampled_data_histogram.svg)](sampled_data_histogram.pdf)

# An equidistant grid is used.
# Adjusting the grid size for sampling can impact the sampling process,
# check `generate` method for details.




# ## Likelihood and Model Fitting

# The concept of unbinned likelihood is central to many statistical models.
# In `AlgebraPDF.jl`, the negative log-likelihood is implemented.
# Creating a negative log-likelihood function for our model and data.

nll = NegativeLogLikelihood(model, data)

# Evaluating the negative log-likelihood gives us an idea of how well the model fits the data.
# The negative log likelihood function depends only on the model parameters, not on the variable `x`.

nll(1.0), nll(0.0), nll(())

# The extended version adds a Poisson factor for the total number of events with expectation given by
# the norm of the model, constraining the model normalization to the number of size of the data sample.

ext = Extended(nll)

# To fit the model to the data, we need to find the parameter values that minimize the extended NLL.
# We start by setting initial values for the parameters, most importantly
# the normalization that is going to be constrained to size of the data set.

starting_values = let
    Nd = length(data)
    default_values = pars(ext)
    @unpack N1, N2 = default_values
    Nsum = N1+N2
    default_values + (N1 = N1/Nsum*Nd, N2 = N2/Nsum*Nd)
end

# When the optimization to fit the model to the data,
# using a reasonable guess for the inverse hessian matrix helps the initial steps of the optimization.
# The diagonal elements of the invesse hessian reflect the typical variation of the parameters.
# The parameter uncertainties in the minimum are often taken as square root of the diagonal elements. 

@time fit = let
    initial_invH = Diagonal([0.001,0.01,0.01,100,100]) .+ eps()
    
    optimize(x->ext(1.1, x), starting_values |> collect,
	    BFGS(; initial_invH = x -> initial_invH,))
end

# Once we have the best-fit parameters, we can update our model and compare it to the original data.

best_model = updatepars(model, NamedTuple{keys(pars(model))}(fit.minimizer))

# Plotting the original data, and the best-fit model helps us evaluate the fitting process.

let
    bins = range(lims(model)..., 100)
    Nd = length(data)

    stephist(data; bins)
    plot!(model, scaletobinneddata(Nd, bins), lab="original")
    plot!(best_model, scaletobinneddata(bins), lab="fit")

    plot!(best_model[2], scaletobinneddata(bins), fill=0, lab="signal")
    plot!(best_model[1], scaletobinneddata(bins), ls=:dash, lab="background")
end
#jl savefig("model_fitting.pdf")
#md savefig("model_fitting.pdf")
#md savefig("model_fitting.svg"); nothing # hide
#md # [![Model Fitting](model_fitting.svg)](model_fitting.pdf)
