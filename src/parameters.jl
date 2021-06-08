const ∅ = NamedTuple()

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
releasepars(ps::NamedTuple, s, v = 0.0) = complain_about_Pars()

#              _|                                      _|      
#    _|_|_|  _|_|_|_|  _|  _|_|  _|    _|    _|_|_|  _|_|_|_|  
#  _|_|        _|      _|_|      _|    _|  _|          _|      
#      _|_|    _|      _|        _|    _|  _|          _|      
#  _|_|_|        _|_|  _|          _|_|_|    _|_|_|      _|_|  
                                                             
struct TwoNamedTuples{R,S}
    free::NamedTuple{R}
    fixed::NamedTuple{S}
end

TwoNamedTuples(t::NamedTuple) = TwoNamedTuples(t,∅)
TwoNamedTuples(ps::TwoNamedTuples) = TwoNamedTuples(freepars(ps),fixedpars(ps))
TwoNamedTuples(; kw...) = TwoNamedTuples((;kw...))

nt(s::Symbol, v = 0.0) = NamedTuple{(s,)}([v])

+(ps::TwoNamedTuples, t2::NamedTuple) = TwoNamedTuples(freepars(ps)+t2, fixedpars(ps))
+(t2::NamedTuple, ps::TwoNamedTuples) = +(ps, t2)
+(ps1::TwoNamedTuples, ps2::TwoNamedTuples) = TwoNamedTuples(freepars(ps1)+freepars(ps2), fixedpars(ps1)+fixedpars(ps2))

# 
freepars(ps::TwoNamedTuples) = getfield(ps, :free)
fixedpars(ps::TwoNamedTuples) = getfield(ps, :fixed)
allpars(ps::TwoNamedTuples) = freepars(ps) + fixedpars(ps)
# 
import Base: getproperty
getproperty(ps::TwoNamedTuples, s::Symbol) = hasproperty(freepars(ps), s) ? getproperty(freepars(ps), s) : getproperty(fixedpars(ps), s)

fixpars(ps::TwoNamedTuples, from_p) = TwoNamedTuples(freepars(ps) - from_p, fixedpars(ps) + from_p)
fixpar(ps::TwoNamedTuples, s::Symbol, v::T where T <: Real = getproperty(freepars(ps),s)) = fixpars(ps, nt(s,v))
releasepars(ps::TwoNamedTuples, sv) =
    TwoNamedTuples(freepars(ps) + sv, fixedpars(ps) - sv)
releasepar(ps::TwoNamedTuples, s::Symbol) =
    TwoNamedTuples(freepars(ps) + nt(s, getproperty(fixedpars(ps),s)), fixedpars(ps) - s)
#
#
updatepars(ps::TwoNamedTuples, from_p::NamedTuple) = TwoNamedTuples(updatepars(freepars(ps), from_p), updatepars(fixedpars(ps), from_p))

selectpars(  p, symb) = NamedTuple{Tuple(symb)}(getproperty.(Ref(p), symb))
subtractpars(p, symb) = Base.structdiff(p, selectpars(p, symb))

copy(pars::TwoNamedTuples, p::NamedTuple) = error("copy(pars, p) works with {pars} ≤ {p}")
copy(pars::NamedTuple, p::TwoNamedTuples) = selectintersect(freepars(pars), allpars(p))
function copy(pars::TwoNamedTuples, pnew::TwoNamedTuples)
    pnewstructure = TwoNamedTuples(
        selectintersect(freepars(pnew), allpars(pars)),
        selectintersect(fixedpars(pnew), allpars(pars)))
    updatepars(pnewstructure, allpars(pnew))
end

