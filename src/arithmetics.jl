
# operation
*(c, d::pdf) = *(d::pdf, c) # commutation
*(d::pdf, c::NamedTuple) = pdf((x;p=∅)->func(d,x;p=p) * getproperty(p, keys(c)[1]);
        p = pars(d) + c, lims = lims(d))
*(d::pdf, c::Number) = pdf((x;p=∅)->func(d,x;p=p) * c;
        p = pars(d), lims = lims(d))
*(d1::pdf, d2::pdf) = pdf((x;p=∅)->func(d1,x;p=p) .* func(d2,x;p=p);
        p = pars(d1) + pars(d2), lims = lims(d1))
#
/(d1::pdf, d2::pdf) = pdf((x;p=∅)->func(d1,x;p=p) ./ func(d2,x;p=p);
    p = pars(d1) + pars(d2), lims = lims(d1))
# 
+(c::NamedTuple, d::pdf) = pdf((x;p=∅)->func(d,x;p=p) + getproperty(p,keys(c)[1]);
        p = pars(d) + c, lims = lims(d))
# pdf + pdf
+(d1::pdf, d2::pdf) = pdf((x;p=∅)->func(d1,x;p=p) + func(d2,x;p=p);
        p = pars(d1)+pars(d2), lims = lims(d1))
#