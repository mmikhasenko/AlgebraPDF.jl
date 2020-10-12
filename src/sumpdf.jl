

randfractionsymbol() = Symbol("fsum_"*randstring(3))

function sumpdf(f1,f2,lims)
    n1 = quadgk(f1,lims...)[1]
    n2 = quadgk(f2,lims...)[1]
    fracsym = randfractionsymbol()
    d = pdf((x;p)->getproperty(p,fracsym)*f1(x)/n1 + (1-getproperty(p,fracsym))*f2(x)/n2;
        lims=lims, p0= NamedTuple{(fracsym,)}(0.5))
    return d
end

function sumpdf(d1::pdf, d2::pdf, sfS::Symbol=randfractionsymbol())
    d1.lims != d2.lims && error("lims are different: ", d1.lims, " != ", d2.lims)
    n1 = integral(d1)
    n2 = integral(d2)
    d = pdf(
        (x;p)->getproperty(p, sfS) .* d1.f(x;p=p) ./ n1 + (1-getproperty(p, sfS)) .* d2.f(x;p=p) ./ n2;
        lims=d1.lims, p0=merge(d1.p0, d2.p0, NamedTuple{(sfS,)}(0.5)))
    return d
end

function sumpdf(d1::pdf, d2::pdf, f1::Float64)
    d1.lims != d2.lims && error("lims are different: ", d1.lims, " != ", d2.lims)
    n1 = integral(d1)
    n2 = integral(d2)
    d = pdf(
        (x;p)->f1 .* d1.f(x;p=p) ./ n1 + (1-f1) .* d2.f(x;p=p) ./ n2;
        lims=d1.lims, p0=merge(d1.p0, d2.p0))
    return d
end

