
###################################################################### 

struct Abs2Func{T<:AbstractFunctionWithParameters} <: AbstractFunctionWithParameters
    f::T
end
func(d::Abs2Func, x::NumberOrTuple; p=freepars(d)) = abs2(func(d.f,x;p))
pars(d::Abs2Func, isfree::Bool) = pars(d.f, isfree)
updatevalueorflag( d::Abs2Func, s::Symbol, isfree::Bool, v=getproperty(pars(d),s)) =
    Abs2Func(updatevalueorflag(d.f,s,isfree,v))
# 
abs2(f::AbstractFunctionWithParameters) = Abs2Func(f)

###################################################################### 

struct LogFunc{T<:AbstractFunctionWithParameters} <: AbstractFunctionWithParameters
    f::T
end
func(d::LogFunc, x::NumberOrTuple; p=freepars(d)) = log(func(d.f,x;p))
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
func(d::ProdFunc, x::NumberOrTuple; p=freepars(d)) = func(d.f1,x;p) * func(d.f2,x;p)
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
func(d::NegativeLogLikelihood, x::NumberOrTuple; p=freepars(d)) = -sum((v>0) ? log(v) : d.nagativepenatly for v in d.f(d.data;p))
pars(d::NegativeLogLikelihood, isfree::Bool) = pars(d.f, isfree)
updatevalueorflag( d::NegativeLogLikelihood, s::Symbol, isfree::Bool, v=getproperty(pars(d),s)) =
    NegativeLogLikelihood(updatevalueorflag(d.f,s,isfree,v))
#
NegativeLogLikelihood(d, data::AbstractArray) = NegativeLogLikelihood(d, data, -1e4)
minussum(d::LogFunc, data::AbstractArray) = NegativeLogLikelihood(d.f, data)


###################################################################### 

struct SumFunc{T<:AbstractFunctionWithParameters,N,V}  <: AbstractFunctionWithParameters
    fs::StaticVector{N,T}
    αs::V
end
function SumFunc(fs::Union{Tuple,AbstractVector}, αs)
    N = length(fs)
    SumFunc(SVector{N}(collect(fs)), αs)
end
# 
length(d::SumFunc{T,N}) where {T,N} = N
pars(d::SumFunc, isfree::Bool) = sum(pars.(d.fs, isfree)) + pars(d.αs, isfree)
getindex(d::SumFunc, i::Number) = d.fs[i]

const SumOfFunc = AlgebraPDF.SumFunc{T} where T<:AbstractFunctionWithParameters
const SumOfPDF = AlgebraPDF.SumFunc{T} where T<:AbstractPDF

function func(d::SumOfFunc, x::NumberOrTuple; p=freepars(d))
    allp = p+fixedpars(d)
    αs_vals = (getproperty(allp,s) for s in keys(d.αs))
    f_vals = func.(d.fs, Ref(x);p)
    return sum(αs_vals .* f_vals)   
end

# 
function func_norm(d::SumOfPDF, x; p=freepars(d)) # suppose to work also for all x <: AbstractVector
    allp = p+fixedpars(d)
    αs_vals = (getproperty(allp,s) for s in keys(d.αs))
    f_vals = [f(x;p) for f in d.fs]
    return sum(αs_vals .* f_vals)    
end
func(d::SumOfPDF, x::NumberOrTuple; p=freepars(d)) = func_norm(d,x;p)
func(d::SumOfPDF, x::AbstractArray; p=freepars(d)) = func_norm(d,x;p)
func(d::SumOfPDF, x::AbstractRange; p=freepars(d)) = func_norm(d,x;p)
lims(d::SumOfPDF) = lims(d.fs[1])

function updatevalueorflag( d::SumFunc, s::Symbol, isfree::Bool, v=getproperty(pars(d),s))
    fs = [ispar(f,s) ? updatevalueorflag(f,s,isfree,v) : f for f in d.fs]
    αs = ispar(d.αs,s) ? updatevalueorflag(d.αs,s,isfree,v) : d.αs
    SumFunc(fs, αs)
end
#

+(f1::AbstractFunctionWithParameters,
  f2::AbstractFunctionWithParameters, αs=(α1=1.0, α2=1.0)) = SumFunc([f1,f2], αs)

-(f1::AbstractFunctionWithParameters,
  f2::AbstractFunctionWithParameters, αs=(α1=1.0, α2=-1.0)) = SumFunc([f1,f2], αs)

###################################################################### 
