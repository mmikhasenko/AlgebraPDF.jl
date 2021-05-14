

randfractionsymbol() = Symbol("fsum_"*randstring(3))

function sumpdf(f1,f2,lims)
    n1 = quadgk(f1,lims...)[1]
    n2 = quadgk(f2,lims...)[1]
    fracsym = randfractionsymbol()
    d = pdf((x;p)->getproperty(p,fracsym)*f1(x)/n1 + (1-getproperty(p,fracsym))*f2(x)/n2;
        lims=lims, p=TwoNamedTuples(nt(fracsym,0.5)))
    return d
end

function sumpdf(d1::pdf, d2::pdf, sfS::Symbol=randfractionsymbol())
    lims(d1) != lims(d2) && error("lims are different: ", lims(d1), " != ", lims(d2))
    n1 = normalizationintegral(d1)
    n2 = normalizationintegral(d2)
    d = pdf(
        (x;p)->getproperty(p, sfS) .* func(d1,x;p) ./ n1 + (1-getproperty(p, sfS)) .* func(d2,x;p) ./ n2;
        lims=lims(d1), p=TwoNamedTuples(pars(d1) + pars(d2) + nt(sfS,0.5)))
    return d
end

function sumpdf(d1::pdf, d2::pdf, f1::Float64)
    lims(d1) != lims(d2) && error("lims are different: ", lims(d1), " != ", lims(d2))
    n1 = normalizationintegral(d1)
    n2 = normalizationintegral(d2)
    d = pdf(
        (x;p)->f1 .* func(d1,x;p) ./ n1 + (1-f1) .* func(d2,x;p) ./ n2;
        lims=lims(d1), p=TwoNamedTuples(pars(d1) + pars(d2)))
    return d
end

