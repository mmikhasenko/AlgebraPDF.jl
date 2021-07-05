

#    _|_|    _|                    _|                                      _|      _|_|_|    _|_|_|    _|_|_|_|  
#  _|    _|  _|_|_|      _|_|_|  _|_|_|_|  _|  _|_|    _|_|_|    _|_|_|  _|_|_|_|  _|    _|  _|    _|  _|        
#  _|_|_|_|  _|    _|  _|_|        _|      _|_|      _|    _|  _|          _|      _|_|_|    _|    _|  _|_|_|    
#  _|    _|  _|    _|      _|_|    _|      _|        _|    _|  _|          _|      _|        _|    _|  _|        
#  _|    _|  _|_|_|    _|_|_|        _|_|  _|          _|_|_|    _|_|_|      _|_|  _|        _|_|_|    _|        
                                                                                                               

abstract type AbstractPDF{N} <: AbstractFunctionWithParameters end  # N is the dimension

# calls
function (d::AbstractPDF)(x; p=freepars(d))
    allp = p+fixedpars(d)
    normalization = normalizationintegral(d; p=allp)
    if normalization ≈ 0.0
        println("Error: normalization ≈ 0!")
        normalization = 1.0
    end
    return func(d,x; p=allp) / normalization
end

# call on 
(d::AbstractPDF)(x, v::AbstractVector) = d(x; p=v2p(v,d))

# 1 dims
normalizationintegral(d::AbstractPDF{1}; p=freepars(d)) =
    quadgk(x->func(d, x; p=p), lims(d)...)[1]
#
# 2 dims
# function normalizationintegral(d::AbstractPDF{2}; p=freepars(d))
#     xmap = x->mapx_to_unit(x,lims(d))
#     curhe((x,f)->f[1]=func(d, xmap(x); p=p), lims(d)...)[1]
# end


"""
noparsnormf(d::AbstractPDF; p=pars(d))

Returns a single-argument lambda-function with parameters fixed to `p` and normalization computed.
""" 
noparsnormf(d::AbstractPDF; p=pars(d)) = (ns=normalizationintegral(d;p=p); (x;kw...)->func(d,x;p=p)/ns)


# by default, I assume that the fields "lims" is present
lims(d::AbstractPDF) = getfield(d, :lims)
updatevalueorflag(d::AbstractPDF, s::Symbol, isfree::Bool, v=getproperty(pars(d),s)) =
    typeof(d)(updatevalueorflag(d.p, s, isfree, v), d.lims)

# other methods
function integral(d::AbstractPDF{1}, lims; p=freepars(d))
    allpars = p+fixedpars(d)
    quadgk(x->func(d,x; p=allpars), lims...)[1] / normalizationintegral(d; p=allpars)
end

#################################################################### 

struct Normalized{T<:AbstractFunctionWithParameters,L} <: AbstractPDF{1}
    lineshape::T
    lims::L
end
lineshape(d::Normalized) = getfield(d, :lineshape)

# two methods to be defined
import Base: getproperty
getproperty(d::Normalized, sym::Symbol) = sym==:p ? pars(lineshape(d)) : getfield(d, sym)
func(d::Normalized, x::NumberOrTuple; p=pars(d)) = func(lineshape(d), x; p)
pars(d::Normalized, isfree::Bool) = pars(lineshape(d), isfree)
updatevalueorflag(d::Normalized, s::Symbol, isfree::Bool, v=getproperty(pars(d),s)) =
    Normalized(updatevalueorflag(lineshape(d), s, isfree, v), d.lims)


# short cuts
# 1 argument
Normalized(f;p,lims) = Normalized(FunctionWithParameters(f; p), lims)

###################################################################### 

fixedshapepdf(f, lims) = Normalized((x;p)->f(x); lims=lims, p=∅)

###################################################################### 

#    _|                                                    _|      _|_|  
#  _|_|_|_|  _|    _|  _|_|_|      _|_|    _|_|_|      _|_|_|    _|      
#    _|      _|    _|  _|    _|  _|_|_|_|  _|    _|  _|    _|  _|_|_|_|  
#    _|      _|    _|  _|    _|  _|        _|    _|  _|    _|    _|      
#      _|_|    _|_|_|  _|_|_|      _|_|_|  _|_|_|      _|_|_|    _|      
#                  _|  _|                  _|                            
#              _|_|    _|                  _|                            


"""
    @makepdftype MyPDF(x;p) = unnormdensity(x, p.a, p.b)

    Expected form of the expression is `f(x;p)` on the left
"""
macro makepdftype(ex)
    # 
    fpx = ex.args[1].args
    name = fpx[1]
    p = fpx[2].args[1]
    (p != :p) && error("expected format f(x;p) = ... " )
    x = fpx[3]
    body = ex.args[2]
    # 
    quote
        struct $name{T,N} <: AbstractPDF{1}
            p::T
            lims::N
        end
        $(esc(name))(;p,lims) = $(esc(name))(p, lims)

        import AlgebraPDF: func, pars, updatevalueorflag
        #
        $(esc(:func))(d::$(esc(name)), $(esc(x))::NumberOrTuple; p=$(esc(:pars))(d)) = $(esc(body))
    end
end
