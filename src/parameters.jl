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

# 
free(ps::Parameters) = collectpars(ps.free)
fixed(ps::Parameters) = collectpars(ps.fixed)
constrained(ps::Parameters) = collectpars(ps.constrained)
# 
collectpars(ps::Parameters) = free(ps) + fixed(ps)

fixpar(ps::Parameters, s::Symbol, v::T where T <: Real = getproperty(ps.free,s)) =
    Parameters(ps.free - s, ps.fixed + nt(s,v), ps.constrained)
releasepar(ps::Parameters, s::Symbol, v::T where T <: Real = getproperty(ps.fixed,s)) =
    Parameters(ps.free + nt(s, v), ps.fixed - s, ps.constrained)
constrainpar(ps::Parameters, s::Symbol, v::T where T <: Real, e::T where T <: Real) =
    Parameters(ps.free, ps.fixed, ps.constrained + nt(s,(v,e)))
unconstrainpar(ps::Parameters, s::Symbol) =
    Parameters(ps.free, ps.fixed, ps.constrained - s)
# 
updatepars(ps::Parameters, from_p) = Parameters(updatepars(free(ps), from_p), updatepars(fixed(ps), from_p), constrained(ps))