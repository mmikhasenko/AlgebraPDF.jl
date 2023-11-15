"""
    subtract(t1::NamedTuple, t2::NamedTuple)
    subtract(t1::NamedTuple, s::Symbol)
    subtract(t1::NamedTuple, ss::AbstractVector{Symbol})
    subtract(t1::NamedTuple, ss::Tuple{Vararg{Symbol}}) 

The function is used to remove the parameters from the NamedTuple that are given by the second argument.
It can be a single symbol, a vector of symbols or a tuple of symbols.

# Examples
```jldoctest
julia> (a=1,b=2,d=1,c=2) - (d=1,)
(a=1,b=2,c=2)

julia> (a=1,b=2,d=1,c=2) - :d
(a=1,b=2,c=2)

julia> (a=1,b=2,d=1,c=2) - [:d, :a]
(b=2,c=2)

julia> (a=1,b=2,d=1,c=2) - (:d, :a)
(b=2,c=2)
```
"""
subtract(t1::NamedTuple, t2::NamedTuple) = Base.structdiff(t1, t2)
subtract(t1::NamedTuple, s::Symbol) = Base.structdiff(t1, nt(s))
subtract(t1::NamedTuple, ss::AbstractVector{Symbol}) = Base.structdiff(t1, NamedTuple{Tuple(ss)}(zeros(length(ss))))
subtract(t1::NamedTuple, ss::Tuple{Vararg{Symbol}}) = Base.structdiff(t1, NamedTuple{ss}(zeros(length(ss))))


#    _|                          _|            
#  _|_|_|_|  _|    _|  _|_|_|    _|    _|_|    
#    _|      _|    _|  _|    _|  _|  _|_|_|_|  
#    _|      _|    _|  _|    _|  _|  _|        
#      _|_|    _|_|_|  _|_|_|    _|    _|_|_|  
#                      _|                      
#                      _|                      


pars(ps::NamedTuple, isfree::Bool) = isfree == true ? ps : ∅
# 
selectintersect(p::NamedTuple, from_p::NamedTuple) = selectpars(from_p, intersect(keys(p), keys(from_p)))

updatevalueorflag(p::NamedTuple, s::Symbol, isfree::Bool, v=getproperty(pars(p), s)) =
    isfree == true ? NamedTuple{keys(p)}(merge(p, nt(s, v))) : # join and replace the old value
    throw(ArgumentError("Not able to fix, release parameters since parameters are held by a NamedTuple!"))

#              _|                                      _|      
#    _|_|_|  _|_|_|_|  _|  _|_|  _|    _|    _|_|_|  _|_|_|_|  
#  _|_|        _|      _|_|      _|    _|  _|          _|      
#      _|_|    _|      _|        _|    _|  _|          _|      
#  _|_|_|        _|_|  _|          _|_|_|    _|_|_|      _|_|  

"""
    FlaggedNamedTuple(t::NamedTuple) = FlaggedNamedTuple(t, ())
    FlaggedNamedTuple(ps::FlaggedNamedTuple) = FlaggedNamedTuple(allpars(ps), whichfixed(ps))
    FlaggedNamedTuple(; kw...) = FlaggedNamedTuple((; kw...))
    
immutable type that holds a NamedTuple of parameters and a tuple of symbols that indicates which parameters are fixed.
The short alias for the `FlaggedNamedTuple` is `Ext`.

# Examples
```jldoctest
julia> ps = FlaggedNamedTuple((a=1,b=2,c=3), (:a, :b))  # a and b are fixed 
FlaggedNamedTuple{(:a, :b)}((a = 1, b = 2, c = 3), (:a, :b))

julia> freepars(ps)
(c = 3,)

julia> fixedpars(ps)
(a = 1, b = 2)

julia> Ext((a = 2.2,)).a
2.2
```
"""
struct FlaggedNamedTuple{R}
    allpars::NamedTuple{R}
    whichfixed::Tuple{Vararg{Symbol}}
end

FlaggedNamedTuple(t::NamedTuple) = FlaggedNamedTuple(t, ())
FlaggedNamedTuple(ps::FlaggedNamedTuple) = FlaggedNamedTuple(allpars(ps), whichfixed(ps))
FlaggedNamedTuple(; kw...) = FlaggedNamedTuple((; kw...))

allpars(ps::FlaggedNamedTuple) = getfield(ps, :allpars)
whichfixed(ps::FlaggedNamedTuple) = getfield(ps, :whichfixed)

keys(ps::FlaggedNamedTuple) = keys(allpars(ps))
keys(ps::FlaggedNamedTuple, isfree::Bool) = isfree == false ? whichfixed(ps) : Base.diff_names(keys(allpars(ps)), whichfixed(ps))
pars(ps::FlaggedNamedTuple, isfree::Bool) = NamedTuple{keys(ps, isfree)}(allpars(ps))

getproperty(ps::FlaggedNamedTuple, s::Symbol) = hasproperty(freepars(ps), s) ? getproperty(freepars(ps), s) : getproperty(fixedpars(ps), s)


# inner working
updatednamedtuple(p::NamedTuple, s::Symbol, v) = NamedTuple{keys(p)}(merge(p, nt(s, v)))
mustinclude(whichfixed::Tuple{Vararg{Symbol}}, s) = s ∈ whichfixed ? whichfixed : (whichfixed..., s)
mustexclude(whichfixed::Tuple{Vararg{Symbol}}, s) = s ∈ whichfixed ? Base.diff_names(whichfixed, (s,)) : whichfixed
function updateflag(whichfixed::Tuple{Vararg{Symbol}}, s::Symbol, isfree::Bool)
    isfree ? mustexclude(whichfixed, s) : mustinclude(whichfixed, s)
end

"""
    updatevalueorflag(p::FlaggedNamedTuple, s::Symbol, isfree::Bool, v=getproperty(pars(p),s))

Implementation of the main update method for `FlaggedNamedTuple` parameters.
"""
function updatevalueorflag(p::FlaggedNamedTuple, s::Symbol, isfree::Bool, v=getproperty(pars(p), s))
    return FlaggedNamedTuple(
        updatednamedtuple(allpars(p), s, v),
        updateflag(whichfixed(p), s, isfree))
end
merge(s1::FlaggedNamedTuple, s2::FlaggedNamedTuple) =
    FlaggedNamedTuple(merge(allpars(s1), allpars(s2)), (whichfixed(s1)..., whichfixed(s2)...))
merge(s1::FlaggedNamedTuple, s2::NamedTuple) = merge(s1, FlaggedNamedTuple(s2))
merge(s1::NamedTuple, s2::FlaggedNamedTuple) = merge(FlaggedNamedTuple(s1), s2)
#
const ParTypes = Union{NamedTuple,FlaggedNamedTuple}
# 
pars(ps::ParTypes) = NamedTuple{keys(ps)}(merge(pars(ps, true), pars(ps, false)))
freepars(d::ParTypes) = pars(d, true)
fixedpars(d::ParTypes) = pars(d, false)
#
