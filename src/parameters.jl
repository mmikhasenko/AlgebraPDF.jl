const ∅ = NamedTuple()

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
                                                             
struct TwoNamedTuples{R}
    allpars::NamedTuple{R}
    whichfixed::Tuple{Vararg{Symbol}}
end
allpars(ps::TwoNamedTuples) = getfield(ps,:allpars)
whichfixed(ps::TwoNamedTuples) = getfield(ps,:whichfixed)

keys(ps::TwoNamedTuples) = keys(allpars(ps))
keys(ps::TwoNamedTuples, isfree::Bool) = isfree==false ? whichfixed(ps) : Base.diff_names(keys(allpars(ps)), whichfixed(ps))
pars(ps::TwoNamedTuples, isfree::Bool) = NamedTuple{keys(ps, isfree)}(allpars(ps))

getproperty(ps::TwoNamedTuples, s::Symbol) = hasproperty(freepars(ps), s) ? getproperty(freepars(ps), s) : getproperty(fixedpars(ps), s)


# inner working
updatednamedtuple(p::NamedTuple, s::Symbol, v) = NamedTuple{keys(p)}(p+nt(s,v))
mustinclude(whichfixed::Tuple{Vararg{Symbol}},s) = s ∈ whichfixed ? whichfixed : (whichfixed...,s)
mustexclude(whichfixed::Tuple{Vararg{Symbol}},s) = s ∈ whichfixed ? Base.diff_names(whichfixed, (s,)) : whichfixed
function updateflag(whichfixed::Tuple{Vararg{Symbol}}, s::Symbol, isfree::Bool)
    isfree ? mustexclude(whichfixed,s) : mustinclude(whichfixed,s)
end

"""
    updatevalueorflag(p::TwoNamedTuples, s::Symbol, isfree::Bool, v=getproperty(pars(p),s))

Implementation of the main update method for `TwoNamedTuples` parameters.
"""
function updatevalueorflag(p::TwoNamedTuples, s::Symbol, isfree::Bool, v=getproperty(pars(p),s))
    return TwoNamedTuples(
        updatednamedtuple(allpars(p),s,v),
        updateflag(whichfixed(p), s, isfree))
end

TwoNamedTuples(t::NamedTuple) = TwoNamedTuples(t,())
TwoNamedTuples(ps::TwoNamedTuples) = TwoNamedTuples(allpars(ps),whichfixed(ps))
TwoNamedTuples(; kw...) = TwoNamedTuples((;kw...))

const ParTypes = Union{NamedTuple,TwoNamedTuples}
#
pars(ps::ParTypes) = NamedTuple{keys(ps)}(pars(ps, true) + pars(ps, false))
freepars(d::ParTypes) = pars(d, true)
fixedpars(d::ParTypes) = pars(d, false)
