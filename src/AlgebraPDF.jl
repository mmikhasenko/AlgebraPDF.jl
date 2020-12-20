module AlgebraPDF
#
using Parameters
using QuadGK
using Optim
using Random
using RecipesBase
using StaticArrays
using Interpolations
using ForwardDiff
using LinearAlgebra
using Measurements

#
import Base: +, *, /

import Optim: minimizer, minimum
import ForwardDiff: hessian

export AdvancedFunction
export ∅
export pdf, npars, p2v, v2p
export integral, integrals
export fixpars, fixpar
export fixedshapepdf
export noparsf, noparsnormf
export collectpars, func, lims
export subtractpars, selectpars, updatepars
include("pdf.jl")

export fit_llh, llh
export fit_llh_with_constraints
export hessian, covmat, invH, invexacthessian
export errors, measurements
export minimizer, minimum
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

export scaletobinneddata
include("plottingrecipe.jl")

export sWeights
include("sWeights.jl")

export xProductPDF
include("multidim.jl")

export aGauss, aBreitWigner, aExp
export aDoubleGaussFixedRatio
export aBreitWignerConvGauss
export aTabulated
include("densities.jl")

export inrange
include("utils.jl")

end # module
