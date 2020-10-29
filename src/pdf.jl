@with_kw struct pdf
    f::Function
    lims::Tuple{Real,Real}
    xdim::Int = 1
    p0::NamedTuple
end
#
# consructors
pdf(f,p0,lims) = pdf(;f=f,lims=lims,p0=p0)
pdf(f;p0,lims) = pdf(;f=f,lims=lims,p0=p0)
#
collectpars(d::pdf) = d.p0

integral(d::pdf; p=d.p0) = quadgk(x->d.f(x; p=p), d.lims...)[1]

# calls
function (d::pdf)(x; p=d.p0, norm_according_to=d)
    normalization = integral(norm_according_to; p=p)
    normalization ≈ 0.0 && error("norm = 0 with p = $(p)!")
    return d.f(x; p=p) / normalization
end
(d::pdf)(x, v) = d(x; p=v2p(v,d))

# properties
npars(d) = length(collectpars(d))
v2p(v,d) = NamedTuple{keys(collectpars(d))}(v)
p2v(p,d) = [getproperty(p, k) for k in keys(collectpars(d))]
p2v(  d) = p2v(collectpars(d), d)


∅ = NamedTuple()
# operation
*(c, d::pdf) = *(d::pdf, c) # commutation
*(d::pdf, c::NamedTuple) = pdf((x;p=∅)->d.f(x;p=p) * getproperty(p, keys(c)[1]);
        p0 = merge(c, d.p0), lims = d.lims)
*(d::pdf, c::Number) = pdf((x;p=∅)->d.f(x;p=p) * c;
        p0 = d.p0, lims = d.lims)
*(d1::pdf, d2::pdf) = pdf((x;p=∅)->d1.f(x;p=p) .* d2.f(x;p=p);
        p0 = merge(d1.p0, d2.p0), lims = d1.lims)
#
/(d1::pdf, d2::pdf) = pdf((x;p=∅)->d1.f(x;p=p) ./ d2.f(x;p=p);
    p0 = merge(d1.p0, d2.p0), lims = d1.lims)
# 
+(c::NamedTuple, d::pdf) = pdf((x;p=∅)->d.f(x;p=p) + getproperty(p,keys(c)[1]);
        p0 = merge(c, d.p0), lims = d.lims)
# pdf + pdf
+(d1::pdf, d2::pdf) = pdf((x;p=∅)->d1.f(x;p=p) + d2.f(x;p=p);
        p0 = merge(d1.p0, d2.p0), lims = d1.lims)
#
# fix parameters
fix_parameters(d::pdf, values::NamedTuple) = pdf((e;p=∅)->d.f(e;p=merge(p,values)); p0=Base.structdiff(d.p0, values), lims=d.lims)
fix_parameters(d::pdf, symb::T where T<:Union{Tuple,Array{Symbol}}) = fix_parameters(d, NamedTuple{Tuple(symb)}(getproperty.(Ref(d.p0),symb)))
#
fixedshapepdf(f, lims) = pdf((x;p=∅)->f(x); lims=lims, p0=∅)

#
noparsf(d::pdf; p=d.p0) = x->d.f(x;p=p)
noparsnormf(d::pdf; p=d.p0) = (ns=integral(d;p=p); x->d.f(x;p=p)/ns)
