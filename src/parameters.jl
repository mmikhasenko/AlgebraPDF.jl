const ∅ = NamedTuple()

+(t1::NamedTuple, t2::NamedTuple) = merge(t1,t2)
-(t1::NamedTuple, t2::NamedTuple) = Base.structdiff(t1,t2)
-(t1::NamedTuple, s::Symbol) = Base.structdiff(t1, nt(s))
-(t1::NamedTuple, ss::Vector{Symbol}) = Base.structdiff(t1, sum(nt.(ss)))
-(t1::NamedTuple, ss::Tuple) = Base.structdiff(t1, sum(nt.(ss)))

nt(s::Symbol, v = 0.0) = NamedTuple{(s,)}([v])

#    _|                          _|            
#  _|_|_|_|  _|    _|  _|_|_|    _|    _|_|    
#    _|      _|    _|  _|    _|  _|  _|_|_|_|  
#    _|      _|    _|  _|    _|  _|  _|        
#      _|_|    _|_|_|  _|_|_|    _|    _|_|_|  
#                      _|                      
#                      _|                      


pars(ps::NamedTuple, isfree) = isfree==true ? ps : ∅
# 
selectintersect(p::NamedTuple, from_p::NamedTuple) = selectpars(from_p, intersect(keys(p), keys(from_p)))

updatevalueorflag(p::NamedTuple, s::Symbol, isfree::Bool, v=getproperty(pars(d),s)) =
    isfree==true ? NamedTuple{keys(p)}(p+nt(s,v)) : # join and replace the old value
    throw(DomainError("Not able to fix, release parameters since parameters are held by a NamedTuple!"))

#              _|                                      _|      
#    _|_|_|  _|_|_|_|  _|  _|_|  _|    _|    _|_|_|  _|_|_|_|  
#  _|_|        _|      _|_|      _|    _|  _|          _|      
#      _|_|    _|      _|        _|    _|  _|          _|      
#  _|_|_|        _|_|  _|          _|_|_|    _|_|_|      _|_|  
                                                             
struct TwoNamedTuples{R}
    allpars::NamedTuple{R}
    whichfixed::Vector{Symbol}
end
allpars(ps::TwoNamedTuples) = getfield(ps,:allpars)
whichfixed(ps::TwoNamedTuples) = getfield(ps,:whichfixed)

keys(ps::TwoNamedTuples, isfree::Bool) = isfree==false ? Tuple(whichfixed(ps)) : Base.diff_names(keys(allpars(ps)), Tuple(whichfixed(ps)))
pars(ps::TwoNamedTuples, isfree::Bool) = NamedTuple{keys(ps, isfree)}(allpars(ps))

getproperty(ps::TwoNamedTuples, s::Symbol) = hasproperty(freepars(ps), s) ? getproperty(freepars(ps), s) : getproperty(fixedpars(ps), s)


# +(ps::TwoNamedTuples, t2::NamedTuple) = TwoNamedTuples(freepars(ps)+t2, fixedpars(ps))
# +(t2::NamedTuple, ps::TwoNamedTuples) = +(ps, t2)
# +(ps1::TwoNamedTuples, ps2::TwoNamedTuples) = TwoNamedTuples(freepars(ps1)+freepars(ps2), fixedpars(ps1)+fixedpars(ps2))
#
updatednamedtuple(p::NamedTuple, s::Symbol, v) = NamedTuple{keys(p)}(p+nt(s,v))

mustinclude(whichfixed::Vector,s) = s ∈ whichfixed ? whichfixed : push!(whichfixed,s)
mustexclude(whichfixed::Vector,s) = s ∈ whichfixed ? filter!(x->x!=s,whichfixed) : whichfixed
# 
function updateflag!(whichfixed::Vector{Symbol}, s::Symbol, isfree::Bool)
    isfree ? mustexclude(whichfixed,s) : mustinclude(whichfixed,s)
end

function updatevalueorflag(p::TwoNamedTuples, s::Symbol, isfree::Bool, v=getproperty(pars(p),s))
    # newp = (v == getproperty(pars(p),s)) ? p : 
    newp = TwoNamedTuples(updatednamedtuple(allpars(p),s,v), copy(whichfixed(p)))
    updateflag!(getfield(newp,:whichfixed), s, isfree)
    return newp
end

TwoNamedTuples(t::NamedTuple) = TwoNamedTuples(t,Symbol[])
TwoNamedTuples(ps::TwoNamedTuples) = TwoNamedTuples(freepars(ps),fixedpars(ps))
TwoNamedTuples(; kw...) = TwoNamedTuples((;kw...))

==(ps1::TwoNamedTuples, ps2::TwoNamedTuples) =
    (freepars(ps1)==freepars(ps2))&&(fixedpars(ps1)==fixedpars(ps2))

const ParTypes = Union{NamedTuple,TwoNamedTuples}
#
pars(ps::ParTypes) = pars(ps, true) + pars(ps, false)
freepars(d::ParTypes) = pars(d, true)
fixedpars(d::ParTypes) = pars(d, false)
