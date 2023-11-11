var documenterSearchIndex = {"docs":
[{"location":"api/#API-Documentation","page":"API","title":"API Documentation","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"The API for the package ","category":"page"},{"location":"api/#Standard-functions","page":"API","title":"Standard functions","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"There is a range of predefined functions  implemented as a struct that holds tuple of parameters","category":"page"},{"location":"api/#Functions","page":"API","title":"Functions","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"AlgebraPDF.noparsnormf\nAlgebraPDF.nt","category":"page"},{"location":"api/#AlgebraPDF.noparsnormf","page":"API","title":"AlgebraPDF.noparsnormf","text":"noparsnormf(d::AbstractPDF; p=pars(d))\n\nReturns a single-argument lambda-function with parameters fixed to p and normalization computed.\n\n\n\n\n\n","category":"function"},{"location":"api/#AlgebraPDF.nt","page":"API","title":"AlgebraPDF.nt","text":"nt(s::Symbol, v = 0.0)\n\nCreates a named tuple (s=v,) where s is a provided symbol, and v is the value. \n\n\n\n\n\n","category":"function"},{"location":"api/#Macros","page":"API","title":"Macros","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"AlgebraPDF.@makepdftype\nAlgebraPDF.@makefuntype","category":"page"},{"location":"api/#AlgebraPDF.@makepdftype","page":"API","title":"AlgebraPDF.@makepdftype","text":"@makepdftype MyPDF(x;p) = unnormdensity(x, p.a, p.b)\n\nExpected form of the expression is f(x;p) on the left\n\n\n\n\n\n","category":"macro"},{"location":"api/#AlgebraPDF.@makefuntype","page":"API","title":"AlgebraPDF.@makefuntype","text":"@makefuntype MyPDF(x;p) = unnormdensity(x, p.a, p.b)\n\nExpected form of the expression is `f(x;p)` on the left\n\n\n\n\n\n","category":"macro"},{"location":"api/#Other-Functions","page":"API","title":"Other Functions","text":"","category":"section"},{"location":"api/","page":"API","title":"API","text":"AlgebraPDF.updatevalueorflag","category":"page"},{"location":"api/#AlgebraPDF.updatevalueorflag","page":"API","title":"AlgebraPDF.updatevalueorflag","text":"updatevalueorflag(p::FlaggedNamedTuple, s::Symbol, isfree::Bool, v=getproperty(pars(p),s))\n\nImplementation of the main update method for FlaggedNamedTuple parameters.\n\n\n\n\n\n","category":"function"},{"location":"#AlgebraPDF.jl","page":"Home","title":"AlgebraPDF.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"AlgebraPDF.jl is a tool to construct and operate with complex parametric probability density functions (PDFs). ","category":"page"},{"location":"","page":"Home","title":"Home","text":"Basic functionality:","category":"page"},{"location":"","page":"Home","title":"Home","text":"Attach default values of parameters to a function\nUpdate, fix, release parameters\nconstructing a complex model object from set of function:\nalgebra of functions with parameters, e.g. f₁ + f₂, abs2(f), or log(f).\nOn-fly normalization\nconstruction of mixed models in the form f₁ PDF₁ + f₂ PDF₂ + f₃ PDF₃.\nconstruction of likelihood function and extended likelihood function\nplotting recipes","category":"page"},{"location":"","page":"Home","title":"Home","text":"note: Note\nAlgebraPDF.jl is experimental, feel free to share experience working with the package. Submit an issue with suggestions.","category":"page"},{"location":"#Table-of-contents","page":"Home","title":"Table of contents","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Pages = [\n    \"tutorial.md\",\n    \"api.md\",\n    \"developing.md\",\n    \"license.md\",\n]\nDepth = 1","category":"page"},{"location":"#Citing-AlgebraPDF.jl","page":"Home","title":"Citing AlgebraPDF.jl","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"When using AlgebraPDF.jl for research, teaching or similar,  please drop a link for your research to the issues. It would give me encouragement to create a citable link on Zenodo, write a JOSS paper. I would also be happy to mention your research in the Documentation website.","category":"page"},{"location":"#Learning-Julia","page":"Home","title":"Learning Julia","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"The Julia website provides many links to introductory videos and written tutorials, e.g. \"Intro to Julia\", Think Julia: How to Think Like a Computer Scientist and \"The Fast Track to Julia\". If you are familiar with MATLAB or Python, you may also want to take a look at the \"MATLAB–Python–Julia cheatsheet\".","category":"page"},{"location":"","page":"Home","title":"Home","text":"The in-depth article Why Numba and Cython are not substitutes for Julia explains how Julia addresses several fundamental challenges inherent to scientific high-performance computing.","category":"page"},{"location":"#Installation","page":"Home","title":"Installation","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"using Pkg\nPkg.add(\"AlgebraPDF\")","category":"page"},{"location":"#Usage","page":"Home","title":"Usage","text":"","category":"section"},{"location":"","page":"Home","title":"Home","text":"Provide examples of how to use your package.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"EditURL = \"tutorial_lit.jl\"","category":"page"},{"location":"tutorial/#Tutorial","page":"Tutorial","title":"Tutorial","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"This tutorial guides the user through the basic example of creating a density function which is a sum of a gaussian signal peak, and exponential background, sampling from the distribution, and optimizating its parameters with an unbinned fit.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"using AlgebraPDF, AlgebraPDF.Parameters\nusing LinearAlgebra, Optim\n\nusing Random\nRandom.seed!(100)\n\nusing Plots\ntheme(:wong, frame=:box, xlab=\"x\", lab=\"\", minorticks=true,\n\tguidefontvalign=:top, guidefonthalign=:right,\n\txlim=(:auto,:auto), ylim=(0,:auto), grid=false)","category":"page"},{"location":"tutorial/#Function-With-Parameters","page":"Tutorial","title":"Function With Parameters","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"The package provides a standard wrapper of a function with parameters f(x;p), where x is function variable, p is a structure that holds parameters with their names. The simplest and the most common case is where x would be a number, and p is a named tuple. For example,","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"myf(x; p=(a=1.1,b=2.2)) = x*p.a + p.b/x ;\nnothing #hide","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"The modele introduces a type FunctionWithParameters, which intends to behave like the myf from the user prospective.","category":"page"},{"location":"tutorial/#A-Gaussian-function","page":"Tutorial","title":"A Gaussian function","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"The gaussian function is constructed calling a specific type FGauss, and giving a tuple of parameters with their default values","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"gaussian = FGauss((μ=1.1, σ=0.9))","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"It is on of the predefined examples of functions with parameters, FGauss <: AbstractFunctionWithParameters. The object is callable like a regular function","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"gaussian(1.1)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"when array is passed, the function is broadcasted","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"gaussian(-1.8:0.9:1)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Default values of the parameters can be accessed with pars, and freepars method.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"pars(gaussian)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"One can provide the key argument p with the named tuple of parameters. These tuple is always used instread of the default values.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"gaussian(0.0; p=(; μ = 1.1, σ = 0.9))","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"The parameters can be adjusted","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"gaussian(0.0; p=(; μ = 0.0, σ = 1.9))","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Similar to the regular fuction, the object can be plotted.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"plot(gaussian, -4, 7, fill=0, α=0.8)\nsavefig(\"gaussian.pdf\")\nsavefig(\"gaussian.svg\"); nothing # hide","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"(Image: gaussian)","category":"page"},{"location":"tutorial/#Normalization","page":"Tutorial","title":"Normalization","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"To turn an arbitrary function to the probability density, one need to introduce normalization. This is done by attaching the range (support) to the function.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"nGaussian = Normalized(gaussian, (-4, 7))","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"The same is archived with the pipeline.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"@assert nGaussian == gaussian |> Normalized((-4, 7))","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"The normalized object has the call method regular function. The normalization is computed on fly. It is constly for a single-value call.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"nGaussian(1.1)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"When calling on iterable collection, the normalization is computed once.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"nGaussian(-1.8:0.9:1)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"As before, the parameters can be updates by passing a named tuple","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"nGaussian(0.0; p=(; μ = 1.1, σ = 0.9))","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"For plotting of the normalized function, one does not need to specify the range.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"plot(nGaussian, fill=0, α=0.7)\nsavefig(\"nGaussian.pdf\")\nsavefig(\"nGaussian.svg\"); nothing # hide","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"(Image: nGaussian)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"This section introduces the exponential function in AlgebraPDF.jl and demonstrates how to normalize it. Understanding how to work with different types of functions and their normalization is crucial in creating complex distributions.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Exponential function is another lineshape defined in the package. We are using the expnential distribution with a slope α to define the background to our gaussian signal.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"exponential = FExp((; α=-0.2))","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Next, we normalize this function over a specific range in the same way we did with the Gaussian function in the previous sections. Normalization is essential to ensure that the function behaves as a probability density function and we can countably measure the background fraction","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"nExponent = exponential |> Normalized((-4, 7))","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Let's plot the normalized exponential function to visualize it.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"plot(nExponent)\nsavefig(\"nExponential.pdf\")\nsavefig(\"nExponential.svg\"); nothing # hide","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"(Image: nExponential)","category":"page"},{"location":"tutorial/#Summation-of-Functions","page":"Tutorial","title":"Summation of Functions","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Now, we can explore how to sum different types of functions or distributions. Summing functions is essential part of the package functionality that allows one to model more complex distributions.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Here, we create a model that is a sum of the previously defined normalized exponential function and the normalized Gaussian function.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"model = FSum([nExponent, nGaussian], (N1=0.85, N2=0.15))","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"The sum of functions, is an object of type FSum, that hold a list of functions and their weights in static vectors of equal sizes. There is an alternative way to formulate an equivalent model","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"@assert model == nExponent * (N1=0.85, ) + nGaussian * (N2=0.15,)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"where the product of the function with a named tuple return FSum object. The summation between two FSum objects is defined. Complementary, individual components with their weight can be accessed by indexing the FSum object as in the following protting code.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"begin\n    plot(model)\n    plot!(model[1], ls=:dash, lab=\"background\")\n    plot!(model[2], fill=0, lab=\"signal\")\nend\nsavefig(\"model_components.pdf\")\nsavefig(\"model_components.svg\"); nothing # hide","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"(Image: Model and Components)","category":"page"},{"location":"tutorial/#Sampling-from-the-model","page":"Tutorial","title":"Sampling from the model","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Sampling from a model is a key operation in statistical modeling. AlgebraPDF.jl implements sampling through the numerical Inversion Method.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Here, we sample data points from our model.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"@time data = AlgebraPDF.rand(model, 10_000)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Plotting the sampled data as a histogram gives us a visual representation of the distribution.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"stephist(data, bins=100)\nsavefig(\"sampled_data_histogram.pdf\")\nsavefig(\"sampled_data_histogram.svg\"); nothing # hide","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"(Image: Sampled Data Histogram)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"An equidistant grid is used. Adjusting the grid size for sampling can impact the sampling process, check generate method for details.","category":"page"},{"location":"tutorial/#Likelihood-and-Model-Fitting","page":"Tutorial","title":"Likelihood and Model Fitting","text":"","category":"section"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"The concept of unbinned likelihood is central to many statistical models. In AlgebraPDF.jl, the negative log-likelihood is implemented. Creating a negative log-likelihood function for our model and data.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"nll = NegativeLogLikelihood(model, data)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Evaluating the negative log-likelihood gives us an idea of how well the model fits the data. The negative log likelihood function depends only on the model parameters, not on the variable x.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"nll(1.0), nll(0.0), nll(())","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"The extended version adds a Poisson factor for the total number of events with expectation given by the norm of the model, constraining the model normalization to the number of size of the data sample.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"ext = Extended(nll)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"To fit the model to the data, we need to find the parameter values that minimize the extended NLL. We start by setting initial values for the parameters, most importantly the normalization that is going to be constrained to size of the data set.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"starting_values = let\n    Nd = length(data)\n    default_values = pars(ext)\n    @unpack N1, N2 = default_values\n    Nsum = N1+N2\n    merge(default_values, (N1 = N1/Nsum*Nd, N2 = N2/Nsum*Nd))\nend","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"When the optimization to fit the model to the data, using a reasonable guess for the inverse hessian matrix helps the initial steps of the optimization. The diagonal elements of the invesse hessian reflect the typical variation of the parameters. The parameter uncertainties in the minimum are often taken as square root of the diagonal elements.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"@time fit = let\n    initial_invH = Diagonal([0.001,0.01,0.01,100,100]) .+ eps()\n\n    optimize(x->ext(1.1, x), starting_values |> collect,\n\t    BFGS(; initial_invH = x -> initial_invH,))\nend","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Once we have the best-fit parameters, we can update our model and compare it to the original data.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"best_model = updatepars(model, NamedTuple{keys(pars(model))}(fit.minimizer))","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"Plotting the original data, and the best-fit model helps us evaluate the fitting process.","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"let\n    bins = range(lims(model)..., 100)\n    Nd = length(data)\n\n    stephist(data; bins)\n    plot!(model, scaletobinneddata(Nd, bins), lab=\"original\")\n    plot!(best_model, scaletobinneddata(bins), lab=\"fit\")\n\n    plot!(best_model[2], scaletobinneddata(bins), fill=0, lab=\"signal\")\n    plot!(best_model[1], scaletobinneddata(bins), ls=:dash, lab=\"background\")\nend\nsavefig(\"model_fitting.pdf\")\nsavefig(\"model_fitting.svg\"); nothing # hide","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"(Image: Model Fitting)","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"","category":"page"},{"location":"tutorial/","page":"Tutorial","title":"Tutorial","text":"This page was generated using Literate.jl.","category":"page"}]
}
