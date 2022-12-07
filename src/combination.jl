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

struct FSum{
        T<:AbstractFunctionWithParameters,N,V}  <: AbstractFunctionWithParameters
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

function dividenorm(d::AbstractFunctionWithParameters, Nsymb::Symbol, lims)
    I = quadgk(x->func(d, x), lims...)[1]
    FSum([d], NamedTuple{(Nsymb,)}(1/I))
end
# lambda version
dividenorm(Nsymb::Symbol, lims) = d->dividenorm(d, Nsymb, lims)

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
