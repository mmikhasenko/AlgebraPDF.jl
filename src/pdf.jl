abstract type AdvancedFunction end

const uTAS = Union{Tuple,Array{Symbol}}
∅ = NamedTuple()
freepars(f::Function) = ∅
# properties
npars(d::T where T<:AdvancedFunction) = length(freepars(d))
v2p(v,d::T where T<:AdvancedFunction) = NamedTuple{keys(freepars(d))}(v)
p2v(p,d::T where T<:AdvancedFunction) = [getproperty(p, k) for k in keys(freepars(d))]
p2v(  d::T where T<:AdvancedFunction) = p2v(freepars(d), d)

fixpars(d::T where T<:AdvancedFunction, s::T where T<:uTAS) =
        fixpars(d, selectpars(freepars(d), s))
fixpar(d::T where T<:AdvancedFunction, s::Symbol) =  fixpars(d, (s,))
fixpar(d::T where T<:AdvancedFunction, s::Symbol, v::T where T<:Real) =  fixpars(d, nt(s,v))

                               
#                  _|      _|_|  
#  _|_|_|      _|_|_|    _|      
#  _|    _|  _|    _|  _|_|_|_|  
#  _|    _|  _|    _|    _|      
#  _|_|_|      _|_|_|    _|      
#  _|                            
#  _|                            

@with_kw struct pdf{T} <: AdvancedFunction
    f::Function
    lims::Tuple{Real,Real}
    xdim::Int = 1
    p::T
end
#
# consructors
pdf(f,p,lims) = pdf(;f=f,lims=lims,p=Parameters(p))
pdf(f;p,lims) = pdf(;f=f,lims=lims,p=Parameters(p))
fixedshapepdf(f, lims) = pdf((x;p=∅)->f(x); lims=lims, p=∅)
#
#
lims(d::pdf) = d.lims
pars(d::pdf) = d.p
freepars(d::pdf) = freepars(d.p)
func(d,x; p=pars(d)) = d.f(x;p=p)
# 
normalizationintegral(d::pdf; p=freepars(d)) = quadgk(x->func(d,x; p=p), lims(d)...)[1]
integral(d::pdf, lims; p=freepars(d)) =
        quadgk(x->func(d,x; p=p), lims...)[1] / normalizationintegral(d; p=p)
#
# calls
function (d::pdf)(x; p=freepars(d), norm_according_to=d)
    allp = pars(d) + p
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
fixpars(d::pdf, s_or_from_p::NamedTuple) = pdf(;f=d.f, lims=d.lims, p=fixpars(pars(d), s_or_from_p))
fixpars(d::pdf{T} where T <: NamedTuple, from_p) = error("fixing parameters on pdf{NamedTuple} is outdated! pdf{Parameters} should be constructed by default.")
#
noparsf(d::pdf; p=pars(d)) = (x;kw...)->func(d,x;p=p)
noparsnormf(d::pdf; p=pars(d)) = (ns=normalizationintegral(d;p=p); (x;kw...)->func(d,x;p=p)/ns)
