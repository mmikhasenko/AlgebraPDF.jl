abstract type FunctionWithParameters end

const uTAS = Union{Tuple,Array{Symbol}}
∅ = NamedTuple()
freepars(f::Function) = ∅
#
freepars(d::T where T<:FunctionWithParameters) = freepars(pars(d))
fixedpars(d::T where T<:FunctionWithParameters) = fixedpars(pars(d))
# 
# properties
npars(d::T where T<:FunctionWithParameters) = length(freepars(d))
v2p(v,d::T where T<:FunctionWithParameters) = NamedTuple{keys(freepars(d))}(v)
p2v(p,d::T where T<:FunctionWithParameters) = [getproperty(p, k) for k in keys(freepars(d))]
p2v(  d::T where T<:FunctionWithParameters) = p2v(freepars(d), d)

fixpars(d::T where T<:FunctionWithParameters, s::T where T<:uTAS) =
        fixpars(d, selectpars(freepars(d), s))
fixpar(d::T where T<:FunctionWithParameters, s::Symbol) =  fixpars(d, (s,))
fixpar(d::T where T<:FunctionWithParameters, s::Symbol, v::T where T<:Real) =  fixpars(d, nt(s,v))

#                  _|      _|_|  
#  _|_|_|      _|_|_|    _|      
#  _|    _|  _|    _|  _|_|_|_|  
#  _|    _|  _|    _|    _|      
#  _|_|_|      _|_|_|    _|      
#  _|                            
#  _|                            

@with_kw struct pdf{T} <: FunctionWithParameters
    f::Function
    lims::Tuple{Real,Real}
    xdim::Int = 1
    p::T
end
#
import Base:copy
copy(d::pdf, p) = pdf(;f=func(d), lims=d.lims, p=p)
# 
# consructors
pdf(f,p,lims) = pdf(;f=f,lims=lims,p=Parameters(p))
pdf(f;p,lims) = pdf(;f=f,lims=lims,p=Parameters(p))
fixedshapepdf(f, lims) = pdf((x;p)->f(x); lims=lims, p=∅)
#
#
lims(d::pdf) = d.lims
pars(d::pdf) = d.p
func(d) = d.f
func(d,x; p) = d.f(x;p=p)
#
# 
normalizationintegral(d::pdf; p=freepars(d)) = quadgk(x->func(d,x; p=p), lims(d)...)[1]
function integral(d::pdf, lims; p=freepars(d))
    allpars = p+fixedpars(d)
    quadgk(x->func(d,x; p=allpars), lims...)[1] / normalizationintegral(d; p=allpars)
end
#
# calls
function (d::pdf)(x; p=freepars(d), norm_according_to=d)
    allp = p+fixedpars(d)
    normalization = normalizationintegral(norm_according_to; p=allp)
    if normalization ≈ 0.0
        println("Error: normalization ≈ 0!")
        normalization = 1.0
    end
#     normalization ≈ 0.0 && error("norm = 0 with p = $(p)!")
    return func(d,x; p=allp) / normalization
end
(d::pdf)(x, v) = d(x; p=v2p(v,d))

# fix parameters
fixpars(d::pdf, args::NamedTuple) = copy(d, fixpars(pars(d), args))
releasepar(d::pdf, args...) = copy(d, releasepar(pars(d), args...))
constrainpar(d::pdf, args...) = copy(d, constrainpar(pars(d), args...))
unconstrainpar(d::pdf, args...) = copy(d, unconstrainpar(pars(d), args...))
#
noparsf(d::pdf; p=pars(d)) = (x;kw...)->func(d,x;p=p)
noparsnormf(d::pdf; p=pars(d)) = (ns=normalizationintegral(d;p=p); (x;kw...)->func(d,x;p=p)/ns)
#
updatepars(d::pdf, from_p::NamedTuple) = pdf(;f=func(d), lims=d.lims, p=updatepars(pars(d), from_p))
#

