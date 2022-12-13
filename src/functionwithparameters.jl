pars(d) = pars(d, true) + pars(d, false)
freepars(d) = pars(d, true)
fixedpars(d) = pars(d, false)
ispar(d, s::Symbol) = (s ∈ keys(pars(d)))
isfreepar(d, s::Symbol) = (s ∈ keys(freepars(d)))
#

abstract type AbstractFunctionWithParameters end
#
nfreepars(d::AbstractFunctionWithParameters) = length(freepars(d))
# 
const ArrayOrRange = Union{AbstractArray,AbstractRange}
func(d::AbstractFunctionWithParameters, x::ArrayOrRange; p=freepars(d)) = func.(Ref(d), x; p=p+fixedpars(p))
# 
(d::AbstractFunctionWithParameters)(x; p=freepars(d)) = func(d,x;p)
(d::AbstractFunctionWithParameters)(x, v::ArrayOrRange) = d(x; p=v2p(v,d))

# 0-dims function (NLL)
(d::AbstractFunctionWithParameters)(; p=freepars(d)) = func(d,();p)

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
# lambda versions
updatepar(s::Symbol, v) = d->updatepar(d, s, v) 
releasepar(s::Symbol) = d->releasepar(d, s) 
fixpar(    s::Symbol) = d->fixpar(d, s) 
fixpar(    s::Symbol, v) = d->fixpar(d, s, v) 


# plural
const SymbolSequenceType = Union{AbstractVector{Symbol}, Tuple{Vararg{Symbol}}}
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
# 
updatepars(sequence::NamedTuple) = d->updatepars(d, sequence)
fixpars(sequence::NamedTuple) = d->fixpars(d, sequence)
# 
# Number <: AbstractFunctionWithParameters
pars(d::Number, isfree::Bool) = ∅
func(d::Number, x::NumberOrTuple; p=∅) = d

# Function <: AbstractFunctionWithParameters
pars(f::Function, isfree::Bool) = ∅
func(f::Function, x::NumberOrTuple; p=∅) = f(x)

# # to be implemented for F <: AbstractFunctionWithParameters
# func(d::F, x::NumberOrTuple; p)
# pars(d::F, isfree::Bool)
# updatevalueorflag( d::F, s::Symbol, isfree::Bool, v=getproperty(pars(d),s))

# default method assumes that the d has d.p, and is a single arg
pars(d::AbstractFunctionWithParameters, isfree::Bool) = pars(d.p, isfree)
updatevalueorflag(d::AbstractFunctionWithParameters, s::Symbol, isfree::Bool, v=getproperty(pars(d),s)) =
    typeof(d)(updatevalueorflag(d.p, s, isfree, v))
#

# wrap julia function
@with_kw struct FunctionWithParameters{T} <: AbstractFunctionWithParameters
    f::Function
    p::T
end
FunctionWithParameters(f;p) = FunctionWithParameters(;f,p)

# two methods to be defined
func(d::FunctionWithParameters, x::NumberOrTuple; p=freepars(d)) = d.f(x; p=p+fixedpars(d))
pars(d::FunctionWithParameters, isfree::Bool) = pars(d.p, isfree)
updatevalueorflag(d::FunctionWithParameters, s::Symbol, isfree::Bool, v=getproperty(pars(d),s)) =
    FunctionWithParameters(;f=d.f, p=updatevalueorflag(d.p,s,isfree,v))
#


"""
    @makefuntype MyPDF(x;p) = unnormdensity(x, p.a, p.b)

    Expected form of the expression is `f(x;p)` on the left
"""
macro makefuntype(ex)
    # 
    fpx = ex.args[1].args
    name = fpx[1]
    p = fpx[2].args[1]
    (p != :p) && error("expected format f(x;p) = ... " )
    x = fpx[3]
    body = ex.args[2]
    #
    quote
        struct $name{T} <: AbstractFunctionWithParameters
            p::T
        end
        $(esc(name))(;p) = $(esc(name))(p)
        # 
        $(esc(:(AlgebraPDF.func)))(d::$(esc(name)), $(esc(x))::NumberOrTuple;
            p=$(esc(:(AlgebraPDF.pars)))(d)) =
                $(esc(body))
    end
end

# parameters conversion
v2p(v,d::AbstractFunctionWithParameters) = NamedTuple{keys(freepars(d))}(v)
p2v(p,d::AbstractFunctionWithParameters) = [getproperty(p, k) for k in keys(freepars(d))]
p2v(  d::AbstractFunctionWithParameters) = p2v(freepars(d), d)

# single-argument lambda-function with fixed parameters
noparsf(d::AbstractFunctionWithParameters; p=pars(d)) = (x;kw...)->func(d,x;p=p)


# numerical normalization
function cumulativefunc(d::AbstractFunctionWithParameters, x1::NumberOrTuple, x2::NumberOrTuple; p=freepars(d))
    allp = p+fixedpars(d)
    return quadgk(x->func(d, x; p=allp), x1, x2)[1]
end
