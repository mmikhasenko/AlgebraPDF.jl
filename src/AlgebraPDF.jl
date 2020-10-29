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

export ∅

export fit_llh, llh
include("fit.jl")

export pdf, npars, p2v, v2p
export integral
export fix_parameters
export fixedshapepdf
export noparsf, noparsnormf
export collectpars
include("pdf.jl")

export MixedModel
export update_namedtuple
export fractions
include("mixedmodel.jl")

export generate
include("generation.jl")

export conv_with_gauss, conv_with_gauss_sampling
include("convolution.jl")

export sumpdf
include("sumpdf.jl")

include("plottingrecipe.jl")

export sWeights
include("sWeights.jl")

export xProductPDF
include("multidim.jl")


end # module
