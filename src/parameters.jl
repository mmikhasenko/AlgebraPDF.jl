struct Parameters{R,S,T}
    free::R
    fixed::S
    constrained::T
end

Pars(; kw...) = Parameters((;kw...))
Parameters(t::NamedTuple) = Parameters(t,∅,∅)
Parameters(ps::Parameters) = Parameters(freepars(ps),fixedpars(ps),constrainedpars(ps))

nt(s::Symbol, v = 0.0) = NamedTuple{(s,)}([v])

import Base:+,-
+(t1::NamedTuple, t2::NamedTuple) = merge(t1,t2)
-(t1::NamedTuple, t2::NamedTuple) = Base.structdiff(t1,t2)
-(t1::NamedTuple, s::Symbol) = Base.structdiff(t1, nt(s))
-(t1::NamedTuple, ss::Vector{Symbol}) = Base.structdiff(t1, sum(nt.(ss)))
-(t1::NamedTuple, ss::Tuple) = Base.structdiff(t1, sum(nt.(ss)))

+(ps::Parameters, t2::NamedTuple) = Parameters(freepars(ps)+t2, fixedpars(ps), constrainedpars(ps))
+(ps1::Parameters, ps2::Parameters) = Parameters(freepars(ps1)+freepars(ps2), fixedpars(ps1)+fixedpars(ps2), constrainedpars(ps1)+constrainedpars(ps2))

# 
freepars(ps::Parameters) = getfield(ps, :free)
fixedpars(ps::Parameters) = getfield(ps, :fixed)
constrainedpars(ps::Parameters) = getfield(ps, :constrained)
# 
import Base: getproperty
getproperty(ps::Parameters, s::Symbol) = hasproperty(freepars(ps), s) ? getproperty(freepars(ps), s) : getproperty(fixedpars(ps), s)

fixpars(ps::Parameters, from_p) = Parameters(freepars(ps) - from_p, fixedpars(ps) + from_p, constrainedpars(ps))
fixpar(ps::Parameters, s::Symbol, v::T where T <: Real = getproperty(freepars(ps),s)) = fixpars(ps, nt(s,v))
# 
releasepar(ps::Parameters, s::Symbol, v::T where T <: Real = getproperty(fixedpars(ps),s)) =
    Parameters(freepars(ps) + nt(s, v), fixedpars(ps) - s, constrainedpars(ps))
constrainpar(ps::Parameters, s::Symbol, v::T where T <: Real, e::T where T <: Real) =
    Parameters(freepars(ps), fixedpars(ps), constrainedpars(ps) + nt(s,(v,e)))
unconstrainpar(ps::Parameters, s::Symbol) =
    Parameters(freepars(ps), fixedpars(ps), constrainedpars(ps) - s)
#
#
selectintersect(p::NamedTuple, from_p::NamedTuple) = selectpars(from_p, intersect(keys(p), keys(from_p)))
updatepars(p::NamedTuple, from_p::NamedTuple) = p + selectintersect(p, from_p)
updatepars(ps::Parameters, from_p::NamedTuple) = Parameters(updatepars(freepars(ps), from_p), updatepars(fixedpars(ps), from_p), constrainedpars(ps))

selectpars(  p, symb) = NamedTuple{Tuple(symb)}(getproperty.(Ref(p), symb))
subtractpars(p, symb) = Base.structdiff(p, selectpars(p, symb))

# NamedTuple
freepars(ps::NamedTuple) = ps
fixedpars(ps::NamedTuple) = ∅
constrainedpars(ps::NamedTuple) = ∅

#
complain_about_Pars() = throw(DomainError("To be able to fix, release and constrain parameters create pdf with the `Pars(a=1.1, b=2.2)` constructor"))
fixpars(ps::NamedTuple, from_p) = complain_about_Pars()
fixpar(ps::NamedTuple, s::Symbol, v::T where T <: Real = 0.0) = complain_about_Pars()
releasepar(ps::NamedTuple, s::Symbol, v::T where T <: Real = 0.0) = complain_about_Pars()
constrainpar(ps::NamedTuple, s::Symbol, v::T where T <: Real, e::T where T <: Real) = complain_about_Pars()
unconstrainpar(ps::NamedTuple, s::Symbol) = complain_about_Pars()