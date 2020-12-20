abstract type AdvancedFunction end

∅ = NamedTuple()
collectpars(f::Function) = ∅
# properties
npars(d::T where T<:AdvancedFunction) = length(collectpars(d))
v2p(v,d::T where T<:AdvancedFunction) = NamedTuple{keys(collectpars(d))}(v)
p2v(p,d::T where T<:AdvancedFunction) = [getproperty(p, k) for k in keys(collectpars(d))]
p2v(  d::T where T<:AdvancedFunction) = p2v(collectpars(d), d)

fixpars(d::T where T<:AdvancedFunction, symb::T where T<:Union{Tuple,Array{Symbol}}) =
        fixpars(d, selectpars(collectpars(d), symb))
fixpar(d::T where T<:AdvancedFunction, symb::Symbol) =  fixpars(d, (symb,))
fixpar(d::T where T<:AdvancedFunction, symb::Symbol, value::Float64) =
        fixpars(d, NamedTuple{(symb,)}(value))
@with_kw struct pdf <: AdvancedFunction
    f::Function
    lims::Tuple{Real,Real}
    xdim::Int = 1
    p::NamedTuple
end
#
# consructors
pdf(f,p,lims) = pdf(;f=f,lims=lims,p=p)
pdf(f;p,lims) = pdf(;f=f,lims=lims,p=p)
fixedshapepdf(f, lims) = pdf((x;p=∅)->f(x); lims=lims, p=∅)
#
#
lims(d::pdf) = d.lims
collectpars(d::pdf) = d.p
func(d,x; p=collectpars(d)) = d.f(x;p=p)
# 
normalizationintegral(d::pdf; p=collectpars(d)) = quadgk(x->func(d,x; p=p), lims(d)...)[1]
integral(d::pdf, lims; p=collectpars(d)) =
        quadgk(x->func(d,x; p=p), lims...)[1] / normalizationintegral(d; p=p)
# 

# calls
function (d::pdf)(x; p=collectpars(d), norm_according_to=d)
    normalization = normalizationintegral(norm_according_to; p=p)
    if normalization ≈ 0.0
        println("Error: normalization ≈ 0!")
        normalization = 1.0
    end
#     normalization ≈ 0.0 && error("norm = 0 with p = $(p)!")
    return func(d,x; p=p) / normalization
end
(d::pdf)(x, v) = d(x; p=v2p(v,d))

# operation
*(c, d::pdf) = *(d::pdf, c) # commutation
*(d::pdf, c::NamedTuple) = pdf((x;p=∅)->func(d,x;p=p) * getproperty(p, keys(c)[1]);
        p = merge(c, collectpars(d)), lims = lims(d))
*(d::pdf, c::Number) = pdf((x;p=∅)->func(d,x;p=p) * c;
        p = collectpars(d), lims = lims(d))
*(d1::pdf, d2::pdf) = pdf((x;p=∅)->func(d1,x;p=p) .* func(d2,x;p=p);
        p = merge(collectpars(d1), collectpars(d2)), lims = lims(d1))
#
/(d1::pdf, d2::pdf) = pdf((x;p=∅)->func(d1,x;p=p) ./ func(d2,x;p=p);
    p = merge(collectpars(d1), collectpars(d2)), lims = lims(d1))
# 
+(c::NamedTuple, d::pdf) = pdf((x;p=∅)->func(d,x;p=p) + getproperty(p,keys(c)[1]);
        p = merge(c, collectpars(d)), lims = lims(d))
# pdf + pdf
+(d1::pdf, d2::pdf) = pdf((x;p=∅)->func(d1,x;p=p) + func(d2,x;p=p);
        p = merge(collectpars(d1), collectpars(d2)), lims = lims(d1))
#
# fix parameters

subtractpars(p, symb) = Base.structdiff(p, selectpars(p, symb))
selectpars(  p, symb) = NamedTuple{Tuple(symb)}(getproperty.(Ref(p), symb))
selectintersect(p, from_p) = selectpars(from_p, intersect(keys(p), keys(from_p)))
updatepars(  p, from_p) = merge(p, selectintersect(p, from_p))
# NamedTuple{Tuple(keys(p))}(
#         [hasproperty(from_p, s) ? getproperty(from_p, s) : getproperty(p, s) for s in keys(p)])
# 
fixpars(d::pdf, pars::NamedTuple) =
    (pars == ∅) ? d : pdf((e;p=∅)->func(d,e; p=merge(p,pars));
        p=subtractpars(collectpars(d), keys(pars)), lims=lims(d))
#
noparsf(d::pdf; p=collectpars(d)) = (x;kw...)->func(d,x;p=p)
noparsnormf(d::pdf; p=collectpars(d)) = (ns=normalizationintegral(d;p=p); (x;kw...)->func(d,x;p=p)/ns)
