abstract type FunctionWithParameters end

∅ = NamedTuple()
freepars(f::Function) = ∅
#
freepars(d::FunctionWithParameters) = freepars(pars(d))
fixedpars(d::FunctionWithParameters) = fixedpars(pars(d))
# 
# properties
npars(d::FunctionWithParameters) = length(freepars(d))
v2p(v,d::FunctionWithParameters) = NamedTuple{keys(freepars(d))}(v)
p2v(p,d::FunctionWithParameters) = [getproperty(p, k) for k in keys(freepars(d))]
p2v(  d::FunctionWithParameters) = p2v(freepars(d), d)

# 
fixpar(d::FunctionWithParameters, s::Symbol) =  fixpars(d, (s,))
fixpar(d::FunctionWithParameters, s::Symbol, v::Real) =  fixpars(d, nt(s,v))
fixpars(d::FunctionWithParameters, s::Union{Tuple,Array{Symbol}}) =
    fixpars(d, selectpars(freepars(d), s))
#
fixpars(       d::FunctionWithParameters, args...) = copy(d, fixpars(       pars(d), args...))
releasepar(    d::FunctionWithParameters, args...) = copy(d, releasepar(    pars(d), args...))
constrainpar(  d::FunctionWithParameters, args...) = copy(d, constrainpar(  pars(d), args...))
unconstrainpar(d::FunctionWithParameters, args...) = copy(d, unconstrainpar(pars(d), args...))
updatepars(    d::FunctionWithParameters, args...) = copy(d, updatepars(    pars(d), args...))


#    _|_|    _|                    _|                                      _|      _|_|_|    _|_|_|    _|_|_|_|  
#  _|    _|  _|_|_|      _|_|_|  _|_|_|_|  _|  _|_|    _|_|_|    _|_|_|  _|_|_|_|  _|    _|  _|    _|  _|        
#  _|_|_|_|  _|    _|  _|_|        _|      _|_|      _|    _|  _|          _|      _|_|_|    _|    _|  _|_|_|    
#  _|    _|  _|    _|      _|_|    _|      _|        _|    _|  _|          _|      _|        _|    _|  _|        
#  _|    _|  _|_|_|    _|_|_|        _|_|  _|          _|_|_|    _|_|_|      _|_|  _|        _|_|_|    _|        
                                                                                                               

abstract type AbstractPDF <: FunctionWithParameters end

normalizationintegral(d::AbstractPDF; p=freepars(d)) = quadgk(x->func(d,x; p=p), lims(d)...)[1]
function integral(d::AbstractPDF, lims; p=freepars(d))
    allpars = p+fixedpars(d)
    quadgk(x->func(d,x; p=allpars), lims...)[1] / normalizationintegral(d; p=allpars)
end
#
# calls
function (d::AbstractPDF)(x; p=freepars(d), norm_according_to=d)
    allp = p+fixedpars(d)
    normalization = normalizationintegral(norm_according_to; p=allp)
    if normalization ≈ 0.0
        println("Error: normalization ≈ 0!")
        normalization = 1.0
    end
#     normalization ≈ 0.0 && error("norm = 0 with p = $(p)!")
    return func(d,x; p=allp) / normalization
end
(d::AbstractPDF)(x, v) = d(x; p=v2p(v,d))
func(d::AbstractPDF,x::AbstractArray{T,N} where {T,N}; p=pars(d)) = func.(Ref(d), x; p=p)
func(d::AbstractPDF,x::AbstractRange{T} where T; p=pars(d)) = func.(Ref(d), x; p=p)

# assumes that the fields "lims" and "p" are present
lims(d::AbstractPDF) = d.lims
pars(d::AbstractPDF) = d.p

                                                                       
#    _|                                                    _|      _|_|  
#  _|_|_|_|  _|    _|  _|_|_|      _|_|    _|_|_|      _|_|_|    _|      
#    _|      _|    _|  _|    _|  _|_|_|_|  _|    _|  _|    _|  _|_|_|_|  
#    _|      _|    _|  _|    _|  _|        _|    _|  _|    _|    _|      
#      _|_|    _|_|_|  _|_|_|      _|_|_|  _|_|_|      _|_|_|    _|      
#                  _|  _|                  _|                            
#              _|_|    _|                  _|                            


macro typepdf(name)
    quote
        struct $name{T,N} <: AbstractPDF
            p::T
            lims::N
        end
        $(esc(name))(;p,lims) = $(esc(name))(Pars(;p...), lims)

        import AlgebraPDF: func
        # func(d::$name, x::Number; p=pars(d)) = $f
    end
end

#                  _|      _|_|  
#  _|_|_|      _|_|_|    _|      
#  _|    _|  _|    _|  _|_|_|_|  
#  _|    _|  _|    _|    _|      
#  _|_|_|      _|_|_|    _|      
#  _|                            
#  _|                            


@with_kw struct pdf{T} <: AbstractPDF
    f::Function
    p::T
    lims::Tuple{Real,Real}
end
pdf(f;p,lims) = pdf(;f=f,lims=lims,p=Parameters(p))
#
# two methods to be defined
func(d::pdf,x::Number; p=pars(d)) = d.f(x;p=p)
copy(d::pdf, p) = pdf(;f=d.f, lims=d.lims, p=p)
#

# 
#
noparsf(d::pdf; p=pars(d)) = (x;kw...)->func(d,x;p=p)
noparsnormf(d::pdf; p=pars(d)) = (ns=normalizationintegral(d;p=p); (x;kw...)->func(d,x;p=p)/ns)


fixedshapepdf(f, lims) = pdf((x;p)->f(x); lims=lims, p=∅)
