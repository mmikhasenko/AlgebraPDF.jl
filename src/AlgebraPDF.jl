module AlgebraPDF
#
using Parameters
using QuadGK
using Optim
using Random
using RecipesBase
using StaticArrays
#
import Base: +, *, /

export AdvancedFunction
export ∅
export pdf, npars, p2v, v2p
export integral
export fixpars, fixpar
export fixedshapepdf
export noparsf, noparsnormf
export collectpars, func, lims
export subtractpars, selectpars, updatepars
include("pdf.jl")

export fit_llh, llh
include("fit.jl")

export MixedModel
export fractions
include("mixedmodel.jl")

export generate
include("generation.jl")

export conv_with_gauss
export conv_with_gauss_sampling
include("convolution.jl")

export sumpdf
include("sumpdf.jl")

include("plottingrecipe.jl")

export sWeights
include("sWeights.jl")

export xProductPDF
include("multidim.jl")

export aGauss, aBreitWigner, aExp
include("densities.jl")


end # module
