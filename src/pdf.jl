abstract type AbstractFunctionWithParameters end

# 
const ∅ = NamedTuple()
freepars(f::Function) = ∅
#
pars(d::AbstractFunctionWithParameters) = d.p
freepars(d::AbstractFunctionWithParameters) = freepars(pars(d))
fixedpars(d::AbstractFunctionWithParameters) = fixedpars(pars(d))
# 
# properties
npars(d::AbstractFunctionWithParameters) = length(freepars(d))
v2p(v,d::AbstractFunctionWithParameters) = NamedTuple{keys(freepars(d))}(v)
p2v(p,d::AbstractFunctionWithParameters) = [getproperty(p, k) for k in keys(freepars(d))]
p2v(  d::AbstractFunctionWithParameters) = p2v(freepars(d), d)

# 
fixpar(d::AbstractFunctionWithParameters, s::Symbol) =  fixpars(d, (s,))
fixpar(d::AbstractFunctionWithParameters, s::Symbol, v::Real) =  fixpars(d, nt(s,v))
fixpars(d::AbstractFunctionWithParameters, s::Union{Tuple,Array{Symbol}}) =
    fixpars(d, selectpars(freepars(d), s))
#
fixpars(       d::AbstractFunctionWithParameters, args...) = copy(d, fixpars(       pars(d), args...))
releasepar(    d::AbstractFunctionWithParameters, args...) = copy(d, releasepar(    pars(d), args...))
constrainpar(  d::AbstractFunctionWithParameters, args...) = copy(d, constrainpar(  pars(d), args...))
unconstrainpar(d::AbstractFunctionWithParameters, args...) = copy(d, unconstrainpar(pars(d), args...))
updatepars(    d::AbstractFunctionWithParameters, args...) = copy(d, updatepars(    pars(d), args...))


#    _|_|    _|                    _|                                      _|      _|_|_|    _|_|_|    _|_|_|_|  
#  _|    _|  _|_|_|      _|_|_|  _|_|_|_|  _|  _|_|    _|_|_|    _|_|_|  _|_|_|_|  _|    _|  _|    _|  _|        
#  _|_|_|_|  _|    _|  _|_|        _|      _|_|      _|    _|  _|          _|      _|_|_|    _|    _|  _|_|_|    
#  _|    _|  _|    _|      _|_|    _|      _|        _|    _|  _|          _|      _|        _|    _|  _|        
#  _|    _|  _|_|_|    _|_|_|        _|_|  _|          _|_|_|    _|_|_|      _|_|  _|        _|_|_|    _|        
                                                                                                               

abstract type AbstractPDF <: AbstractFunctionWithParameters end

normalizationintegral(d::AbstractPDF; p=freepars(d)) =
    quadgk(x->func(d, x; p=p), lims(d)...)[1]
function integral(d::AbstractPDF, lims; p=freepars(d))
    allpars = p+fixedpars(d)
    quadgk(x->func(d,x; p=allpars), lims...)[1] / normalizationintegral(d; p=allpars)
end
#
# calls
function (d::AbstractPDF)(x; p=freepars(d))
    allp = p+fixedpars(d)
    normalization = normalizationintegral(d; p=allp)
    if normalization ≈ 0.0
        println("Error: normalization ≈ 0!")
        normalization = 1.0
    end
#     normalization ≈ 0.0 && error("norm = 0 with p = $(p)!")
    return func(d,x; p=allp) / normalization
end

(d::AbstractPDF)(x, v) = d(x; p=v2p(v,d))
func(d::AbstractPDF, x::AbstractArray; p=pars(d)) = func.(Ref(d), x; p=p)
func(d::AbstractPDF, x::AbstractRange; p=pars(d)) = func.(Ref(d), x; p=p)

# assumes that the fields "lims" and "p" are present
lims(d::AbstractPDF) = d.lims


#    _|                                                    _|      _|_|  
#  _|_|_|_|  _|    _|  _|_|_|      _|_|    _|_|_|      _|_|_|    _|      
#    _|      _|    _|  _|    _|  _|_|_|_|  _|    _|  _|    _|  _|_|_|_|  
#    _|      _|    _|  _|    _|  _|        _|    _|  _|    _|    _|      
#      _|_|    _|_|_|  _|_|_|      _|_|_|  _|_|_|      _|_|_|    _|      
#                  _|  _|                  _|                            
#              _|_|    _|                  _|                            

"""
    @typepdf MyPDF(x;p) = unnormdensity(x, p.a, p.b)

    Expected form of the expression is `f(x;p)` on the left
"""
macro typepdf(ex)
    # 
    fpx = ex.args[1].args
    name = fpx[1]
    p = fpx[2].args[1]
    (p != :p) && error("expected format f(x;p) = ... " )
    x = fpx[3]
    body = ex.args[2]
    # 
    quote
        struct $name{T,N} <: AbstractPDF
            p::T
            lims::N
        end
        $(esc(name))(;p,lims) = $(esc(name))(Parameters(p), lims)

        import AlgebraPDF: func
        import Base: copy
        # 
        $(esc(:func))(d::$(esc(name)), $(esc(x))::Number; p=$(esc(:pars))(d)) = $(esc(body))
        $(esc(:copy))(d::$(esc(name)), p) = $(esc(name))(;p=p,lims=$(esc(:lims))(d))
    end
end


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
FunctionWithParameters(f;p) = FunctionWithParameters(;f=f,p=Parameters(p))

# two methods to be defined
func(d::FunctionWithParameters, x::Number; p=pars(d)) = d.f(x; p=p)
copy(d::FunctionWithParameters, p) = FunctionWithParameters(;f=d.f, p=p)
#

#################################################################### 

@with_kw struct pdf{T<:AbstractFunctionWithParameters,N} <: AbstractPDF
    lineshape::T
    lims::N
end
pdf(f;p,lims) = pdf(;
    lineshape = FunctionWithParameters(f; p=Parameters(p)),
        lims = lims)

# two methods to be defined
import Base:getproperty
getproperty(d::pdf, sym::Symbol) = sym==:p ? pars(d.lineshape) : getfield(d, sym)
func(d::pdf, x::Number; p=pars(d)) = func(d.lineshape, x; p=p)
copy(d::pdf, p) = pdf(; lineshape=copy(d.lineshape,p), lims=lims(d))
#

# 
#
noparsf(d::pdf; p=pars(d)) = (x;kw...)->func(d,x;p=p)
noparsnormf(d::pdf; p=pars(d)) = (ns=normalizationintegral(d;p=p); (x;kw...)->func(d,x;p=p)/ns)


fixedshapepdf(f, lims) = pdf((x;p)->f(x); lims=lims, p=∅)


# pdf((e;p)->e^2+p.a; lims=(-1,2), p=(a=1.0,))
# # 
# pdf(
#     FunctionWithParameters((e;p)->e^2+p.a; p=(a=1.0,)),
#     (-1,2))