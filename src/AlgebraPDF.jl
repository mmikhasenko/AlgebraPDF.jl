module AlgebraPDF
#
# using Base: NamedTuple
using ForwardDiff: getindex
using Parameters
using QuadGK
# 
using PyCall
using Optim
# 
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

import Base: +,-,*,==
import Base: getproperty
import Base: getindex
import Base: abs2, log
import Base: keys
import Base: length
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

export FSum, FAbs2, FLog
export NegativeLogLikelihood
export Extended
export FSum
include("arithmetics.jl")

const pdf = Normalized
export pdf

export fit
export OptimHesseApprox
export MigradAndHesse
export minussum
export obj2nt
include("fit.jl")
include("minuit_interface.jl")
include("optim_interface.jl")

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
