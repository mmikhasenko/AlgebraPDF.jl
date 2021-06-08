
#

abstract type AbstractFunctionWithParameters end
#
pars(d::AbstractFunctionWithParameters) = d.p  # default most-common behavior
freepars(d::AbstractFunctionWithParameters) = freepars(pars(d))
fixedpars(d::AbstractFunctionWithParameters) = fixedpars(pars(d))
# 
npars(d::AbstractFunctionWithParameters) = length(freepars(d))



# immutable requires copying
fixpars(       d::AbstractFunctionWithParameters, args...) = copy(d, fixpars(       pars(d), args...))
releasepars(   d::AbstractFunctionWithParameters, args...) = copy(d, releasepars(   pars(d), args...))
updatepars(    d::AbstractFunctionWithParameters, args...) = copy(d, updatepars(    pars(d), args...))

# mutable
fixpars!(      d::AbstractFunctionWithParameters, args...) = error("Does not find mutable implementation for this type")
releasepars!(  d::AbstractFunctionWithParameters, args...) = error("Does not find mutable implementation for this type")
updatepars!(   d::AbstractFunctionWithParameters, args...) = error("Does not find mutable implementation for this type")



# Number <: AbstractFunctionWithParameters
pars(d::Number) = ∅
func(d::Number, x::Number; p=∅) = d

# Function <: AbstractFunctionWithParameters
pars(f::Function) = ∅
func(f::Function, x::Number; p=∅) = f(x)



# sugar
fixpar(d::AbstractFunctionWithParameters, s::Symbol) =  fixpars(d, (s,))
fixpar(d::AbstractFunctionWithParameters, s::Symbol, v::Real) =  fixpars(d, nt(s,v))
fixpars(d::AbstractFunctionWithParameters, s::Union{Tuple,Array{Symbol}}) =
    fixpars(d, selectpars(freepars(d), s))
#


"""
    @newfunc MyPDF(x;p) = unnormdensity(x, p.a, p.b)

    Expected form of the expression is `f(x;p)` on the left
"""
macro newfunc(ex)
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

        import AlgebraPDF: func
        import Base: copy
        # 
        $(esc(:func))(d::$(esc(name)), $(esc(x))::Number; p=$(esc(:pars))(d)) = $(esc(body))
        $(esc(:copy))(d::$(esc(name)), p) = $(esc(name))(;p=copy(d.p, p))
    end
end


struct Abs2Func{T<:AbstractFunctionWithParameters} <: AbstractFunctionWithParameters
    f::T
end
func(d::Abs2Func, x::Number; p=pars(d)) = abs2(func(d.f,x;p))
copy(d::Abs2Func, p) = Abs2Func(copy(d.f,p))
pars(d::Abs2Func) = pars(d.f)
# 
import Base: abs2
abs2(f::AbstractFunctionWithParameters) = Abs2Func(f)

struct SumFunc{
        T1<:AbstractFunctionWithParameters,
        T2<:AbstractFunctionWithParameters,
        V}  <: AbstractFunctionWithParameters
    f1::T1
    f2::T2
    α::V
end
func(d::SumFunc, x::Number; p=pars(d)) = func(d.f1,x;p) + getproperty(p,keys(d.α)[1])*func(d.f2,x;p)
copy(d::SumFunc, p) = SumFunc(copy(d.f1, p), copy(d.f2, p), copy(d.α, p))
pars(d::SumFunc) = pars(d.f1) + pars(d.f2) + d.α
# 
+(f1::AbstractFunctionWithParameters,
  f2::AbstractFunctionWithParameters, α = (α=1.0,)) = SumFunc(f1,f2, α)


#                  _|      _|_|  
#  _|_|_|      _|_|_|    _|      
#  _|    _|  _|    _|  _|_|_|_|  
#  _|    _|  _|    _|    _|      
#  _|_|_|      _|_|_|    _|      
#  _|                            
#  _|                            


@with_kw struct FunctionWithParameters{T} <: AbstractFunctionWithParameters
    f::Function
    p::T
end
FunctionWithParameters(f;p) = FunctionWithParameters(;f,p)

# two methods to be defined
func(d::FunctionWithParameters, x::Number; p=pars(d)) = d.f(x; p)
copy(d::FunctionWithParameters, p) = FunctionWithParameters(;f=d.f, p)
#



# parameters conversion
v2p(v,d::AbstractFunctionWithParameters) = NamedTuple{keys(freepars(d))}(v)
p2v(p,d::AbstractFunctionWithParameters) = [getproperty(p, k) for k in keys(freepars(d))]
p2v(  d::AbstractFunctionWithParameters) = p2v(freepars(d), d)
