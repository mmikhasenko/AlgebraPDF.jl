module AlgebraPDF
#
# using Base: NamedTuple
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
const ∅ = NamedTuple()
export ∅
# 
const NumberOrTuple = Union{Number,Tuple{Vararg{Number}}}
export NumberOrTuple

import Base:+,-,==
import Base: getproperty
import Base: abs2, log
import Base: keys
import Optim: minimizer, minimum
import ForwardDiff: hessian
# 
# 
export ParTypes
include("parameters.jl")
const Ext = FlaggedNamedTuple
export Ext, FlaggedNamedTuple

# 
export updatepar, updatepars
export fixpar, fixpars
export releasepar, releasepars
# 
export @makefuntype
export AbstractFunctionWithParameters
export FunctionWithParameters
export SumFunc, Abs2Func, LogFunc
# 
export func
export nfreepars
export freepars, fixedpars
export pars
export p2v, v2p
# 
include("functionwithparameters.jl")

# 
export @makepdftype
export AbstractPDF
export Normalized
export integral, integrals
export fixedshapepdf
export noparsf, noparsnormf
export lims
export normalizationintegral
# 
include("pdf.jl")
# 
const pdf = Normalized
export pdf


export fit_llh, llh
export fit_llh_with_constraints
export hessian, covmat, invH, invexacthessian
export errors, measurements
export minimizer, minimum
include("fit.jl")

export MixedModel
export fractions, fractionvalues
export Ncomp
include("mixedmodel.jl")

export generate
include("generation.jl")

export convGauss
include("convolution.jl")

export scaletobinneddata
include("plottingrecipe.jl")

export FGauss, FBreitWigner, FExp
export FPowExp, FPol
export FDoubleGaussFixedRatio
export FBreitWignerConvGauss
export FTabulated
include("densities.jl")

export nt
export inrange
include("utils.jl")

end # module
