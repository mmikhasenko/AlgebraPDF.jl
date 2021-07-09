
###################################################################### 

struct Abs2Func{T<:AbstractFunctionWithParameters} <: AbstractFunctionWithParameters
    f::T
end
func(d::Abs2Func, x::NumberOrTuple; p=pars(d)) = abs2(func(d.f,x;p))
pars(d::Abs2Func, isfree::Bool) = pars(d.f, isfree)
updatevalueorflag( d::Abs2Func, s::Symbol, isfree::Bool, v=getproperty(pars(d),s)) =
    Abs2Func(updatevalueorflag(d.f,s,isfree,v))
# 
abs2(f::AbstractFunctionWithParameters) = Abs2Func(f)

###################################################################### 

struct LogFunc{T<:AbstractFunctionWithParameters} <: AbstractFunctionWithParameters
    f::T
end
func(d::LogFunc, x::NumberOrTuple; p=pars(d)) = log(func(d.f,x;p))
pars(d::LogFunc, isfree::Bool) = pars(d.f, isfree)
updatevalueorflag( d::LogFunc, s::Symbol, isfree::Bool, v=getproperty(pars(d),s)) =
    LogFunc(updatevalueorflag(d.f,s,isfree,v))
# 
log(f::AbstractFunctionWithParameters) = LogFunc(f)

###################################################################### 

struct ProdFunc{
        T1<:AbstractFunctionWithParameters,
        T2<:AbstractFunctionWithParameters} <: AbstractFunctionWithParameters
    f1::T1
    f2::T2
end
func(d::ProdFunc, x::NumberOrTuple; p=pars(d)) = func(d.f1,x;p) * func(d.f2,x;p)
pars(d::ProdFunc, isfree::Bool) = pars(d.f1, isfree) + pars(d.f2, isfree)
updatevalueorflag( d::ProdFunc, s::Symbol, isfree::Bool, v=getproperty(pars(d),s)) =
    ProdFunc(
        ispar(d.f1,s) ? updatevalueorflag(d.f1,s,isfree,v) : d.f1,
        ispar(d.f2,s) ? updatevalueorflag(d.f2,s,isfree,v) : d.f2)
# 
*(f1::AbstractFunctionWithParameters, f2::AbstractFunctionWithParameters) = ProdFunc(f1,f2)

###################################################################### 

struct NegativeLogLikelihood{T<:AbstractFunctionWithParameters, D<:AbstractArray} <: AbstractFunctionWithParameters
    f::T
    data::D
    nagativepenatly::Float64
end
func(d::NegativeLogLikelihood, x::NumberOrTuple; p=pars(d)) = -sum((v>0) ? log(v) : d.nagativepenatly for v in d.f(d.data;p))
pars(d::NegativeLogLikelihood, isfree::Bool) = pars(d.f, isfree)
updatevalueorflag( d::NegativeLogLikelihood, s::Symbol, isfree::Bool, v=getproperty(pars(d),s)) =
    NegativeLogLikelihood(updatevalueorflag(d.f,s,isfree,v))
#
NegativeLogLikelihood(d, data::AbstractArray) = NegativeLogLikelihood(d, data, -1e4)
minussum(d::LogFunc, data::AbstractArray) = NegativeLogLikelihood(d.f, data)

###################################################################### 

struct SumFunc{
        T1<:AbstractFunctionWithParameters,
        T2<:AbstractFunctionWithParameters,
        V}  <: AbstractFunctionWithParameters
    f1::T1
    f2::T2
    αs::V
end
function func(d::SumFunc, x::NumberOrTuple; p=pars(d))
    α1, α2 = (getproperty(p,s) for s in keys(d.αs))
    α1*func(d.f1,x;p) + α2*func(d.f2,x;p)
end
pars(d::SumFunc, isfree::Bool) = pars(d.f1, isfree) + pars(d.f2, isfree) + pars(d.αs, isfree)
updatevalueorflag( d::SumFunc, s::Symbol, isfree::Bool, v=getproperty(pars(d),s)) =
    SumFunc(
        ispar(d.f1,s) ? updatevalueorflag(d.f1,s,isfree,v) : d.f1,
        ispar(d.f2,s) ? updatevalueorflag(d.f2,s,isfree,v) : d.f2,
        ispar(d.αs,s) ? updatevalueorflag(d.αs,s,isfree,v) : d.αs)
#
+(f1::AbstractFunctionWithParameters,
  f2::AbstractFunctionWithParameters, αs=(α1=1.0, α2=1.0)) = SumFunc(f1,f2, αs)

-(f1::AbstractFunctionWithParameters,
  f2::AbstractFunctionWithParameters, αs=(α1=1.0, α2=-1.0)) = SumFunc(f1,f2, αs)
###################################################################### 