+(t1::NamedTuple, t2::NamedTuple) = merge(t1,t2)
-(t1::NamedTuple, t2::NamedTuple) = Base.structdiff(t1,t2)
-(t1::NamedTuple, s::Symbol) = Base.structdiff(t1, nt(s))
-(t1::NamedTuple, ss::AbstractVector{Symbol}) = Base.structdiff(t1, sum(nt.(ss)))
-(t1::NamedTuple, ss::Tuple{Vararg{Symbol}}) = Base.structdiff(t1, sum(nt.(ss)))

#    _|                          _|            
#  _|_|_|_|  _|    _|  _|_|_|    _|    _|_|    
#    _|      _|    _|  _|    _|  _|  _|_|_|_|  
#    _|      _|    _|  _|    _|  _|  _|        
#      _|_|    _|_|_|  _|_|_|    _|    _|_|_|  
#                      _|                      
#                      _|                      


pars(ps::NamedTuple, isfree::Bool) = isfree==true ? ps : ∅
# 
selectintersect(p::NamedTuple, from_p::NamedTuple) = selectpars(from_p, intersect(keys(p), keys(from_p)))

updatevalueorflag(p::NamedTuple, s::Symbol, isfree::Bool, v=getproperty(pars(p),s)) =
    isfree==true ? NamedTuple{keys(p)}(p+nt(s,v)) : # join and replace the old value
    throw(ArgumentError("Not able to fix, release parameters since parameters are held by a NamedTuple!"))

#              _|                                      _|      
#    _|_|_|  _|_|_|_|  _|  _|_|  _|    _|    _|_|_|  _|_|_|_|  
#  _|_|        _|      _|_|      _|    _|  _|          _|      
#      _|_|    _|      _|        _|    _|  _|          _|      
#  _|_|_|        _|_|  _|          _|_|_|    _|_|_|      _|_|  
                                                             
struct FlaggedNamedTuple{R}
    allpars::NamedTuple{R}
    whichfixed::Tuple{Vararg{Symbol}}
end
allpars(ps::FlaggedNamedTuple) = getfield(ps,:allpars)
whichfixed(ps::FlaggedNamedTuple) = getfield(ps,:whichfixed)

keys(ps::FlaggedNamedTuple) = keys(allpars(ps))
keys(ps::FlaggedNamedTuple, isfree::Bool) = isfree==false ? whichfixed(ps) : Base.diff_names(keys(allpars(ps)), whichfixed(ps))
pars(ps::FlaggedNamedTuple, isfree::Bool) = NamedTuple{keys(ps, isfree)}(allpars(ps))

getproperty(ps::FlaggedNamedTuple, s::Symbol) = hasproperty(freepars(ps), s) ? getproperty(freepars(ps), s) : getproperty(fixedpars(ps), s)


# inner working
updatednamedtuple(p::NamedTuple, s::Symbol, v) = NamedTuple{keys(p)}(p+nt(s,v))
mustinclude(whichfixed::Tuple{Vararg{Symbol}},s) = s ∈ whichfixed ? whichfixed : (whichfixed...,s)
mustexclude(whichfixed::Tuple{Vararg{Symbol}},s) = s ∈ whichfixed ? Base.diff_names(whichfixed, (s,)) : whichfixed
function updateflag(whichfixed::Tuple{Vararg{Symbol}}, s::Symbol, isfree::Bool)
    isfree ? mustexclude(whichfixed,s) : mustinclude(whichfixed,s)
end

"""
    updatevalueorflag(p::FlaggedNamedTuple, s::Symbol, isfree::Bool, v=getproperty(pars(p),s))

Implementation of the main update method for `FlaggedNamedTuple` parameters.
"""
function updatevalueorflag(p::FlaggedNamedTuple, s::Symbol, isfree::Bool, v=getproperty(pars(p),s))
    return FlaggedNamedTuple(
        updatednamedtuple(allpars(p),s,v),
        updateflag(whichfixed(p), s, isfree))
end

FlaggedNamedTuple(t::NamedTuple) = FlaggedNamedTuple(t,())
FlaggedNamedTuple(ps::FlaggedNamedTuple) = FlaggedNamedTuple(allpars(ps),whichfixed(ps))
FlaggedNamedTuple(; kw...) = FlaggedNamedTuple((;kw...))
#
const ParTypes = Union{NamedTuple,FlaggedNamedTuple}
# 
pars(ps::ParTypes) = NamedTuple{keys(ps)}(pars(ps, true) + pars(ps, false))
freepars(d::ParTypes) = pars(d, true)
fixedpars(d::ParTypes) = pars(d, false)
