struct Parameters{R,S,T}
    free::R
    fixed::S
    constrained::T
end

Parameters(t::NamedTuple) = Parameters(t,∅,∅)
Parameters(ps::Parameters) = Parameters(free(ps),fixed(ps),constrained(ps))

nt(s::Symbol, v = 0.0) = NamedTuple{(s,)}([v])

import Base:+,-
+(t1::NamedTuple, t2::NamedTuple) = merge(t1,t2)
-(t1::NamedTuple, t2::NamedTuple) = Base.structdiff(t1,t2)
-(t1::NamedTuple, s::Symbol) = Base.structdiff(t1, nt(s))
-(t1::NamedTuple, ss::Vector{Symbol}) = Base.structdiff(t1, sum(nt.(ss)))
-(t1::NamedTuple, ss::Tuple) = Base.structdiff(t1, sum(nt.(ss)))

+(ps::Parameters, t2::NamedTuple) = Parameters(free(ps)+t2, fixed(ps), constrained(ps))
+(ps1::Parameters, ps2::Parameters) = Parameters(free(ps1)+free(ps2), fixed(ps1)+fixed(ps2), constrained(ps1)+constrained(ps2))

# 
free(ps::Parameters) = getfield(ps, :free)
fixed(ps::Parameters) = getfield(ps, :fixed)
constrained(ps::Parameters) = getfield(ps, :constrained)
# 
freepars(t::NamedTuple) = t
freepars(ps::Parameters) = free(ps)

import Base: getproperty
getproperty(ps::Parameters, s::Symbol) = hasproperty(free(ps), s) ? getproperty(free(ps), s) : getproperty(fixed(ps), s)

fixpars(ps::Parameters, from_p) = Parameters(free(ps) - from_p, fixed(ps) + from_p, constrained(ps))
fixpar(ps::Parameters, s::Symbol, v::T where T <: Real = getproperty(free(ps),s)) = fixpars(ps, nt(s,v))
# 
releasepar(ps::Parameters, s::Symbol, v::T where T <: Real = getproperty(fixed(ps),s)) =
    Parameters(free(ps) + nt(s, v), fixed(ps) - s, constrained(ps))
constrainpar(ps::Parameters, s::Symbol, v::T where T <: Real, e::T where T <: Real) =
    Parameters(free(ps), fixed(ps), constrained(ps) + nt(s,(v,e)))
unconstrainpar(ps::Parameters, s::Symbol) =
    Parameters(free(ps), fixed(ps), constrained(ps) - s)
#
#
selectintersect(p, from_p) = selectpars(from_p, intersect(keys(p), keys(from_p)))
updatepars(  p, from_p) = merge(p, selectintersect(p, from_p))
updatepars(ps::Parameters, from_p) = Parameters(updatepars(free(ps), from_p), updatepars(fixed(ps), from_p), constrained(ps))

selectpars(  p, symb) = NamedTuple{Tuple(symb)}(getproperty.(Ref(p), symb))
subtractpars(p, symb) = Base.structdiff(p, selectpars(p, symb))
