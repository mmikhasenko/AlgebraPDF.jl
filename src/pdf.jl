

#    _|_|    _|                    _|                                      _|      _|_|_|    _|_|_|    _|_|_|_|  
#  _|    _|  _|_|_|      _|_|_|  _|_|_|_|  _|  _|_|    _|_|_|    _|_|_|  _|_|_|_|  _|    _|  _|    _|  _|        
#  _|_|_|_|  _|    _|  _|_|        _|      _|_|      _|    _|  _|          _|      _|_|_|    _|    _|  _|_|_|    
#  _|    _|  _|    _|      _|_|    _|      _|        _|    _|  _|          _|      _|        _|    _|  _|        
#  _|    _|  _|_|_|    _|_|_|        _|_|  _|          _|_|_|    _|_|_|      _|_|  _|        _|_|_|    _|        
                                                                                                               

abstract type AbstractPDF{N} <: AbstractFunctionWithParameters end  # N is the dimension

# calls
function (d::AbstractPDF)(x; p=freepars(d))
    normalization = normalizationintegral(d; p)
    if normalization ≈ 0.0
        println("Error: normalization ≈ 0!")
        normalization = 1.0
    end
    allp = p+fixedpars(d)
    return func(d,x; p=allp) / normalization
end

# 1 dims
normalizationintegral(d::AbstractPDF; p=freepars(d)) =
    cumulativefunc(d, lims(d)...; p)


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
    return cumulativefunc(d, lims...; p) / normalizationintegral(d; p)
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
# 
function normalizationintegral(d::Normalized; p=freepars(d))
    return cumulativefunc(lineshape(d), lims(d)...; p)
end

# short cuts
# 1 argument
Normalized(lims::L) where L = f->Normalized(f, lims)

###################################################################### 

fixedshapepdf(f, lims) = Normalized(FunctionWithParameters((x;p)->f(x); p=∅), lims)

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

        $(esc(:(AlgebraPDF.func)))(d::$(esc(name)), $(esc(x))::NumberOrTuple;
            p=$(esc(:(AlgebraPDF.pars)))(d)) =
                $(esc(body))
    end
end
