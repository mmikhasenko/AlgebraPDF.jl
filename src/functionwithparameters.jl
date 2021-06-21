
pars(d) = pars(d, true) + pars(d, false)
freepars(d) = pars(d, true)
fixedpars(d) = pars(d, false)
ispar(d, s::Symbol) = (s ∈ keys(pars(d)))
isfreepar(d, s::Symbol) = (s ∈ keys(freepars(d)))
#

abstract type AbstractFunctionWithParameters end
#
npars(d::AbstractFunctionWithParameters) = length(freepars(d))
# methods that call `updatevalueorflag`

updateisfree(d::AbstractFunctionWithParameters, s::Symbol, isfree::Bool) =
    ispar(d,s) ? updatevalueorflag(d,s,isfree) : d
updatevalue( d::AbstractFunctionWithParameters, s::Symbol, v) =
    ispar(d,s) ? updatevalueorflag(d, s, isfreepar(d, s), v) : d
#
# singular
updatepar( d::AbstractFunctionWithParameters, s::Symbol, v) = updatevalue(d, s, v)
releasepar(d::AbstractFunctionWithParameters, s::Symbol) = updateisfree(d, s, true )
fixpar(    d::AbstractFunctionWithParameters, s::Symbol) = updateisfree(d, s, false)
fixpar(    d::AbstractFunctionWithParameters, s::Symbol, v) = updateisfree(updatepar(d, s, v), s, false)

# plural
const SymbolSequenceType = Union{Vector{Symbol}, Tuple}
function fixpars(d::AbstractFunctionWithParameters, sequence::SymbolSequenceType)
    dnew = d
    for s in sequence
        dnew = fixpar(dnew, s)
    end
    return dnew
end

function releasepars(d::AbstractFunctionWithParameters, sequence::SymbolSequenceType)
    dnew = d
    for s in sequence
        dnew = releasepar(dnew, s)
    end
    return dnew
end

function updatepars(d::AbstractFunctionWithParameters, sequence::NamedTuple)
    dnew = d
    for (s,v) in zip(keys(sequence), sequence)
        dnew = updatepar(dnew, s, v)
    end
    return dnew
end
fixpars(d::AbstractFunctionWithParameters, sequence::NamedTuple) = 
    fixpars(updatepars(d, sequence), keys(sequence))

# Number <: AbstractFunctionWithParameters
pars(d::Number, isfree::Bool) = ∅
func(d::Number, x::Number; p=∅) = d

# Function <: AbstractFunctionWithParameters
pars(f::Function, isfree::Bool) = ∅
func(f::Function, x::Number; p=∅) = f(x)

# # to be implemented for F <: AbstractFunctionWithParameters
# func(d::F, x::Number; p)
# pars(d::F, isfree::Bool)
# updatevalueorflag( d::F, s::Symbol, isfree::Bool, v=getproperty(pars(d),s))

###################################################################### 

struct Abs2Func{T<:AbstractFunctionWithParameters} <: AbstractFunctionWithParameters
    f::T
end
func(d::Abs2Func, x::Number; p=pars(d)) = abs2(func(d.f,x;p))
pars(d::Abs2Func, isfree) = pars(d.f, isfree)
updatevalueorflag( d::Abs2Func, s::Symbol, isfree::Bool, v=getproperty(pars(d),s)) =
    Abs2Func(updatevalueorflag(f,s,isfree,v))
# 
abs2(f::AbstractFunctionWithParameters) = Abs2Func(f)

###################################################################### 

struct SumFunc{
        T1<:AbstractFunctionWithParameters,
        T2<:AbstractFunctionWithParameters,
        V}  <: AbstractFunctionWithParameters
    f1::T1
    f2::T2
    α2::V
end
func(d::SumFunc, x::Number; p=pars(d)) = func(d.f1,x;p) + getproperty(p,keys(d.α2)[1])*func(d.f2,x;p)
pars(d::SumFunc) = pars(d.f1) + pars(d.f2) + d.α2
updatevalueorflag( d::SumFunc, s::Symbol, isfree::Bool, v=getproperty(pars(d),s)) =
    SumFunc(
        ispar(d.f1,s) ? updatevalueorflag(d.f1,s,isfree,v) : d.f1,
        ispar(d.f2,s) ? updatevalueorflag(d.f2,s,isfree,v) : d.f2,
        ispar(d.α2,s) ? updatevalueorflag(d.α2,s,isfree,v) : d.α2)
#
+(f1::AbstractFunctionWithParameters,
  f2::AbstractFunctionWithParameters, α2 = (α2=1.0,)) = SumFunc(f1,f2, α2)

###################################################################### 
# wrap julia function
@with_kw struct FunctionWithParameters{T} <: AbstractFunctionWithParameters
    f::Function
    p::T
end
FunctionWithParameters(f;p) = FunctionWithParameters(;f,p)
==(d1::FunctionWithParameters, d2::FunctionWithParameters) =
    (d1.f==d2.f)&&(d1.p==d2.p)

# two methods to be defined
func(d::FunctionWithParameters, x::Number; p=pars(d)) = d.f(x; p)
pars(d::FunctionWithParameters, isfree::Bool) = pars(d.p, isfree)
updatevalueorflag(d::FunctionWithParameters, s::Symbol, isfree::Bool, v=getproperty(pars(d),s)) =
    FunctionWithParameters(;f=d.f, p=updatevalueorflag(d.p,s,isfree,v))
#


# """
#     @newfunc MyPDF(x;p) = unnormdensity(x, p.a, p.b)

#     Expected form of the expression is `f(x;p)` on the left
# """
# macro newfunc(ex)
#     # 
#     fpx = ex.args[1].args
#     name = fpx[1]
#     p = fpx[2].args[1]
#     (p != :p) && error("expected format f(x;p) = ... " )
#     x = fpx[3]
#     body = ex.args[2]
#     # 
#     quote
#         struct $name{T} <: AbstractFunctionWithParameters
#             p::T
#         end
#         $(esc(name))(;p) = $(esc(name))(p)

#         import AlgebraPDF: func
#         import Base: copy
#         # 
#         $(esc(:func))(d::$(esc(name)), $(esc(x))::Number; p=$(esc(:pars))(d)) = $(esc(body))
#         $(esc(:copy))(d::$(esc(name)), p) = $(esc(name))(;p=copy(d.p, p))
#     end
# end

# parameters conversion
v2p(v,d::AbstractFunctionWithParameters) = NamedTuple{keys(freepars(d))}(v)
p2v(p,d::AbstractFunctionWithParameters) = [getproperty(p, k) for k in keys(freepars(d))]
p2v(  d::AbstractFunctionWithParameters) = p2v(freepars(d), d)


