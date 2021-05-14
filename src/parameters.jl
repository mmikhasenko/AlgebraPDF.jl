import Base:+,-
+(t1::NamedTuple, t2::NamedTuple) = merge(t1,t2)
-(t1::NamedTuple, t2::NamedTuple) = Base.structdiff(t1,t2)
-(t1::NamedTuple, s::Symbol) = Base.structdiff(t1, nt(s))
-(t1::NamedTuple, ss::Vector{Symbol}) = Base.structdiff(t1, sum(nt.(ss)))
-(t1::NamedTuple, ss::Tuple) = Base.structdiff(t1, sum(nt.(ss)))


#    _|                          _|            
#  _|_|_|_|  _|    _|  _|_|_|    _|    _|_|    
#    _|      _|    _|  _|    _|  _|  _|_|_|_|  
#    _|      _|    _|  _|    _|  _|  _|        
#      _|_|    _|_|_|  _|_|_|    _|    _|_|_|  
#                      _|                      
#                      _|                      

selectintersect(p::NamedTuple, from_p::NamedTuple) = selectpars(from_p, intersect(keys(p), keys(from_p)))

freepars(ps::NamedTuple) = ps
fixedpars(ps::NamedTuple) = ∅
# 
copy(pars::NamedTuple, pnew::NamedTuple) = selectintersect(pars, pnew)

updatepars(p::NamedTuple, from_p::NamedTuple) = p + selectintersect(p, from_p)
# methods that throw the error
complain_about_Pars() = throw(DomainError("Not able to fix, release parameters since parameters are held by a NamedTuple!"))
fixpars(ps::NamedTuple, from_p)        = complain_about_Pars()
fixpar(ps::NamedTuple, s, v = 0.0)     = complain_about_Pars()
releasepar(ps::NamedTuple, s, v = 0.0) = complain_about_Pars()

#              _|                                      _|      
#    _|_|_|  _|_|_|_|  _|  _|_|  _|    _|    _|_|_|  _|_|_|_|  
#  _|_|        _|      _|_|      _|    _|  _|          _|      
#      _|_|    _|      _|        _|    _|  _|          _|      
#  _|_|_|        _|_|  _|          _|_|_|    _|_|_|      _|_|  
                                                             
struct Parameters{R,S}
    free::R
    fixed::S
end

Pars(; kw...) = Parameters((;kw...))
Parameters(t::NamedTuple) = Parameters(t,∅)
Parameters(ps::Parameters) = Parameters(freepars(ps),fixedpars(ps))

nt(s::Symbol, v = 0.0) = NamedTuple{(s,)}([v])

+(ps::Parameters, t2::NamedTuple) = Parameters(freepars(ps)+t2, fixedpars(ps))
+(t2::NamedTuple, ps::Parameters) = +(ps, t2)
+(ps1::Parameters, ps2::Parameters) = Parameters(freepars(ps1)+freepars(ps2), fixedpars(ps1)+fixedpars(ps2))

# 
freepars(ps::Parameters) = getfield(ps, :free)
fixedpars(ps::Parameters) = getfield(ps, :fixed)
allpars(ps::Parameters) = freepars(ps) + fixedpars(ps)
# 
import Base: getproperty
getproperty(ps::Parameters, s::Symbol) = hasproperty(freepars(ps), s) ? getproperty(freepars(ps), s) : getproperty(fixedpars(ps), s)

fixpars(ps::Parameters, from_p) = Parameters(freepars(ps) - from_p, fixedpars(ps) + from_p)
fixpar(ps::Parameters, s::Symbol, v::T where T <: Real = getproperty(freepars(ps),s)) = fixpars(ps, nt(s,v))
releasepar(ps::Parameters, s::Symbol, v::T where T <: Real = getproperty(fixedpars(ps),s)) =
    Parameters(freepars(ps) + nt(s, v), fixedpars(ps) - s)
#
#
updatepars(ps::Parameters, from_p::NamedTuple) = Parameters(updatepars(freepars(ps), from_p), updatepars(fixedpars(ps), from_p))

selectpars(  p, symb) = NamedTuple{Tuple(symb)}(getproperty.(Ref(p), symb))
subtractpars(p, symb) = Base.structdiff(p, selectpars(p, symb))

copy(pars::Parameters, p::NamedTuple) = error("copy(pars, p) works with {pars} ≤ {p}")
copy(pars::NamedTuple, p::Parameters) = selectintersect(freepars(pars), allpars(p))
function copy(pars::Parameters, pnew::Parameters)
    pnewstructure = Parameters(
        selectintersect(freepars(pnew), allpars(pars)),
        selectintersect(fixedpars(pnew), allpars(pars)))
    updatepars(pnewstructure, allpars(pnew))
end

