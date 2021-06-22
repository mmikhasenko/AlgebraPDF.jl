

#    _|_|    _|                    _|                                      _|      _|_|_|    _|_|_|    _|_|_|_|  
#  _|    _|  _|_|_|      _|_|_|  _|_|_|_|  _|  _|_|    _|_|_|    _|_|_|  _|_|_|_|  _|    _|  _|    _|  _|        
#  _|_|_|_|  _|    _|  _|_|        _|      _|_|      _|    _|  _|          _|      _|_|_|    _|    _|  _|_|_|    
#  _|    _|  _|    _|      _|_|    _|      _|        _|    _|  _|          _|      _|        _|    _|  _|        
#  _|    _|  _|_|_|    _|_|_|        _|_|  _|          _|_|_|    _|_|_|      _|_|  _|        _|_|_|    _|        
                                                                                                               

abstract type AbstractPDF <: AbstractFunctionWithParameters end

normalizationintegral(d::AbstractPDF; p=freepars(d)) =
    quadgk(x->func(d, x; p), lims(d)...)[1]
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
        $(esc(name))(;p,lims) = $(esc(name))(p, lims)

        import AlgebraPDF: func, pars, updatevalueorflag
        # 
        $(esc(:func))(d::$(esc(name)), $(esc(x))::Number; p=$(esc(:pars))(d)) = $(esc(body))
        $(esc(:pars))(d::$(esc(name)), isfree::Bool) = $(esc(:pars))(d.p)
        $(esc(:updatevalueorflag))(d::$(esc(name)), s::Symbol, isfree::Bool, v=getproperty($(esc(:pars))(d),s)) = 
            $(esc(name))(;p=$(esc(:updatevalueorflag))(d.p, s, isfree, v), lims=d.lims)
    end
end


#################################################################### 

@with_kw struct pdf{T<:AbstractFunctionWithParameters,N} <: AbstractPDF
    lineshape::T
    lims::N
end
pdf(f;p,lims) = pdf(;
    lineshape = FunctionWithParameters(f; p),
        lims = lims)

# two methods to be defined
import Base:getproperty
getproperty(d::pdf, sym::Symbol) = sym==:p ? pars(d.lineshape) : getfield(d, sym)
func(d::pdf, x::Number; p=pars(d)) = func(d.lineshape, x; p)
pars(d::pdf, isfree::Bool) = pars(d.lineshape, isfree)
updatevalueorflag(d::pdf, s::Symbol, isfree::Bool, v=getproperty(pars(d),s)) =
    pdf(lineshape=updatevalueorflag(d.lineshape, s, isfree, v), lims=d.lims)

#

# 
#
noparsf(d::AbstractPDF; p=pars(d)) = (x;kw...)->func(d,x;p=p)
noparsnormf(d::AbstractPDF; p=pars(d)) = (ns=normalizationintegral(d;p=p); (x;kw...)->func(d,x;p=p)/ns)

fixedshapepdf(f, lims) = pdf((x;p)->f(x); lims=lims, p=∅)
