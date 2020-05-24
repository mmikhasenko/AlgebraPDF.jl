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
# properties
npars(p::pdf) = length(p.p0)
v2p(v,pdf) = NamedTuple{keys(pdf.p0)}(v)
p2v(p,pdf) = [getproperty(p, k) for k in keys(pdf.p0)]
#
# calls
(d::pdf)(x; p=d.p0, norm_according_to=d) = d.f(x; p=p) / quadgk(x->norm_according_to.f(x; p=p), d.lims...)[1]
(d::pdf)(x, v) = d(x; p=v2p(v,d))
# components
(d::pdf)(x, v)

# operation
*(c, d::pdf) = *(d::pdf, c) # commutation
*(d::pdf, c::NamedTuple) = pdf((x;p)->d.f(x;p=p) * getproperty(p, keys(c)[1]);
        p0 = merge(c, d.p0), lims = d.lims)
*(d::pdf, c::Number) = pdf((x;p)->d.f(x;p=p) * c;
        p0 = d.p0, lims = d.lims)
#
+(c::NamedTuple, d::pdf) = pdf((x;p)->d.f(x;p=p) + getproperty(p,keys(c)[1]);
        p0 = merge(c, d.p0), lims = d.lims)
# pdf + pdf
+(d1::pdf, d2::pdf) = pdf((x;p)->d1.f(x;p=p) + d2.f(x;p=p);
        p0 = merge(d1.p0, d2.p0), lims = d1.lims)
#
