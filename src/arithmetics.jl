
###################################################################### 

struct FAbs2{T<:AbstractFunctionWithParameters} <: AbstractFunctionWithParameters
    f::T
end
func(d::FAbs2, x::NumberOrTuple; p=freepars(d)) = abs2(func(d.f,x;p))
pars(d::FAbs2, isfree::Bool) = pars(d.f, isfree)
updatevalueorflag(d::FAbs2, s::Symbol, isfree::Bool, v=getproperty(pars(d),s)) =
    FAbs2(updatevalueorflag(d.f,s,isfree,v))
# 
abs2(f::AbstractFunctionWithParameters) = FAbs2(f)

###################################################################### 

struct FLog{T<:AbstractFunctionWithParameters} <: AbstractFunctionWithParameters
    f::T
end
func(d::FLog, x::NumberOrTuple; p=freepars(d)) = log(func(d.f,x;p))
pars(d::FLog, isfree::Bool) = pars(d.f, isfree)
updatevalueorflag(d::FLog, s::Symbol, isfree::Bool, v=getproperty(pars(d),s)) =
    FLog(updatevalueorflag(d.f,s,isfree,v))
# 
log(f::AbstractFunctionWithParameters) = FLog(f)


###################################################################### 

struct NegativeLogLikelihood{T<:AbstractFunctionWithParameters, D<:AbstractArray} <: AbstractFunctionWithParameters
    f::T
    data::D
    nagativepenatly::Float64
end
func(d::NegativeLogLikelihood, x::NumberOrTuple; p=freepars(d)) = -sum((v>0) ? log(v) : d.nagativepenatly for v in d.f(d.data;p))
pars(d::NegativeLogLikelihood, isfree::Bool) = pars(d.f, isfree)
updatevalueorflag(d::NegativeLogLikelihood, s::Symbol, isfree::Bool, v=getproperty(pars(d),s)) =
    NegativeLogLikelihood(updatevalueorflag(d.f,s,isfree,v), d.data, d.nagativepenatly)
#
NegativeLogLikelihood(d, data::AbstractArray) = NegativeLogLikelihood(d, data, -1e4)
minussum(d::FLog, data::AbstractArray) = NegativeLogLikelihood(d.f, data)


###################################################################### 

struct ChiSq{M,T<:NumberOrTuple,V<:Number} <: AbstractFunctionWithParameters
    f::M
    xv::Vector{T}
    yv::Vector{V}
    dyv::Vector{V}
end
func(d::ChiSq, x::NumberOrTuple; p=pars(d)) =
    sum((d.yv - d.f(d.xv; p)) .^ 2 ./ d.dyv .^ 2)
# 
pars(d::ChiSq, isfree::Bool) = pars(d.f, isfree)
updatevalueorflag(d::ChiSq, s::Symbol, isfree::Bool, v=getproperty(pars(d),s)) =
    ChiSq(updatevalueorflag(d.f,s,isfree,v), d.xv, d.yv)
#
ChiSq(f, xv::Vector{<:NumberOrTuple}, yv::Vector{V}) where V<:Number =
    ChiSq(f, xv, yv, ones(V,length(yv)))
# 
model(χ²::ChiSq) = χ².f

