
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

struct FProd{
        T1<:AbstractFunctionWithParameters,
        T2<:AbstractFunctionWithParameters} <: AbstractFunctionWithParameters
    f1::T1
    f2::T2
end
func(d::FProd, x::NumberOrTuple; p=freepars(d)) = func(d.f1,x;p) * func(d.f2,x;p)
pars(d::FProd, isfree::Bool) = pars(d.f1, isfree) + pars(d.f2, isfree)
updatevalueorflag(d::FProd, s::Symbol, isfree::Bool, v=getproperty(pars(d),s)) =
    FProd(
        ispar(d.f1,s) ? updatevalueorflag(d.f1,s,isfree,v) : d.f1,
        ispar(d.f2,s) ? updatevalueorflag(d.f2,s,isfree,v) : d.f2)
# 
*(f1::AbstractFunctionWithParameters, f2::AbstractFunctionWithParameters) = FProd(f1,f2)

###################################################################### 

struct FSum{T<:AbstractFunctionWithParameters,N,V}  <: AbstractFunctionWithParameters
    fs::StaticVector{N,T}
    αs::V
end
function FSum(fs::Union{Tuple,AbstractVector}, αs)
    N = length(fs)
    FSum(SVector{N}(collect(fs)), αs)
end
# 
length(d::FSum{T,N}) where {T,N} = N
pars(d::FSum, isfree::Bool) = sum(pars.(d.fs, isfree)) + pars(d.αs, isfree)
function getindex(d::FSum, i::Number)
    α_sym = keys(d.αs)[i]
    α_val = getproperty(d.αs, α_sym)
    FSum([d.fs[i]], nt(α_sym, α_val))
end
const FSumFunc = AlgebraPDF.FSum{T} where T<:AbstractFunctionWithParameters
const FSumPDF = AlgebraPDF.FSum{T} where T<:AbstractPDF

function func(d::FSumFunc, x::NumberOrTuple; p=freepars(d))
    allp = p+fixedpars(d)
    α_vals = (getproperty(allp,s) for s in keys(d.αs))
    f_vals = func.(d.fs, Ref(x);p)
    return sum(α_vals .* f_vals)   
end

# 
function func_norm(d::FSumPDF, x; p=freepars(d)) # suppose to work also for all x <: AbstractVector
    allp = p+fixedpars(d)
    α_vals = (getproperty(allp,s) for s in keys(d.αs))
    f_vals = [f(x;p) for f in d.fs]
    return sum(α_vals .* f_vals)
end
func(d::FSumPDF, x::NumberOrTuple; p=freepars(d)) = func_norm(d,x;p)
func(d::FSumPDF, x::AbstractArray; p=freepars(d)) = func_norm(d,x;p)
func(d::FSumPDF, x::AbstractRange; p=freepars(d)) = func_norm(d,x;p)
lims(d::FSumPDF) = lims(d.fs[1])

function updatevalueorflag(d::FSum, s::Symbol, isfree::Bool, v=getproperty(pars(d),s))
    fs = [ispar(f,s) ? updatevalueorflag(f,s,isfree,v) : f for f in d.fs]
    αs = ispar(d.αs,s) ? updatevalueorflag(d.αs,s,isfree,v) : d.αs
    FSum(fs, αs)
end
#
function integral(d::FSumPDF; p=freepars(d.αs))
    allp = p+fixedpars(d.αs)
    sum(getproperty(allp,s) for s in keys(d.αs))
end

function integral(d::FSumPDF, lims; p=freepars(d))
    allp = p+fixedpars(d)
    α_vals = (getproperty(allp,s) for s in keys(d.αs))
    f_vals = integral.(d.fs, Ref(lims); p)
    return sum(α_vals .* f_vals)
end

+(f1::AbstractFunctionWithParameters,
  f2::AbstractFunctionWithParameters; p=(α1=1.0, α2=1.0)) = FSum([f1,f2], p)

-(f1::AbstractFunctionWithParameters,
  f2::AbstractFunctionWithParameters; p=(α1=1.0, α2=-1.0)) = FSum([f1,f2], p)

+(f1::FSum, f2::FSum) = FSum([f1.fs...,f2.fs...], f1.αs+f2.αs)

*(f::AbstractFunctionWithParameters, c::ParTypes) = FSum([f], c)
*(c::ParTypes, f::AbstractFunctionWithParameters) = FSum([f], c)

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

struct Extended{T <: AlgebraPDF.FSumPDF} <: AbstractFunctionWithParameters
    nll::NegativeLogLikelihood{T}
end

pars(d::Extended, isfree::Bool) = pars(d.nll, isfree)
updatevalueorflag(d::Extended, s::Symbol, isfree::Bool, v=getproperty(pars(d),s)) =
    Extended(updatevalueorflag(d.nll,s,isfree,v))
# 
function func(d::Extended, x::NumberOrTuple; p=freepars(d))
    nll = func(d.nll, x; p)
    μ = integral(d.nll.f; p)
    penalty = μ
    return nll + penalty
end

# model property
model(nll::NegativeLogLikelihood) = nll.f
model(enll::Extended) = model(enll.nll)
