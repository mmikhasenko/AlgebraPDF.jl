module AlgebraPDF
#
using Parameters
using QuadGK
#
using Random
using RecipesBase
using StaticArrays
using Interpolations
using LinearAlgebra
using Measurements
using SpecialFunctions

#
const ∅ = NamedTuple()
export ∅
#
const NumberOrTuple = Union{Number,Tuple{Vararg{Number}}}
export NumberOrTuple

import Base: +, -, *, ==, merge
import Base: getproperty
import Base: getindex
import Base: abs2, log
import Base: keys
import Base: length
#
#
export ParTypes
export subtract
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
export integral
export fixedshapepdf
export noparsf, noparsnormf
export lims
export normalizationintegral
#
include("pdf.jl")
#
export FSum, FAbs2, FLog
export dividenorm
export NegativeLogLikelihood, minussum
export Extended
export ChiSq
include("arithmetics.jl")

const pdf = Normalized
export pdf

export generate
include("generation.jl")

export convGauss
include("convolution.jl")

export scaletobinneddata
include("plottingrecipe.jl")

export FGauss, FBreitWigner, FExp
export FPowExp, FPol
export FLeftSideCrystalBall, FRightSideCrystalBall, FDoubleSideCrystalBall
export FDoubleGaussFixedRatio
export FBreitWignerConvGauss
export FTabulated
include("densities.jl")

export nt
export inrange
include("utils.jl")

end # module
