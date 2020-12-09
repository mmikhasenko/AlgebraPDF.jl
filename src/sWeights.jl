
function Wmatrix(dS, dB, f, lims)
    Wsb = [quadgk(x->di(x)*dj(x)/f(x), lims...)[1]
        for (di,dj) in Iterators.product([dS,dB], [dS,dB])]
    return Wsb
end

function Wmatrix(pdfS::pdf, pdfB::pdf, f)
    nS,nB = normalizationintegral(pdfS),normalizationintegral(pdfB)
    dS = noparsnormf(pdfS)
    dB = noparsnormf(pdfB)
    #
    return Wmatrix(dS, dB, f, lims(pdfS))
end

"""
    sWeights(pdfS::pdf, pdfB::pdf, fraction_of_signal<:Real)

Returns a pair of functions of the discriminating variable `(sWeights_snl(x), sWeights_bgd(x))`.

"""
function sWeights(pdfS::pdf, pdfB::pdf, fraction_of_signal::T where T<:Real)
    ds = sumpdf(pdfS, pdfB, :_f1_sW)
    ds0 = fixpars(ds, (_f1_sW=fraction_of_signal,))
    # 
    Wm = Wmatrix(pdfS, pdfB, noparsf(ds0))
    αm = inv(Wm)
    # 
    αS = αm[:,1] ./ abs(sum(αm[:,1]))
    αB = αm[:,2] ./ abs(sum(αm[:,2]))
    # 
    numerator_snl = fixpars(ds, (_f1_sW=αS[1],))
    sWeights_snl = numerator_snl / ds0 * fraction_of_signal
    # 
    numerator_bgd = fixpars(ds, (_f1_sW=αB[1],))
    sWeights_bgd = numerator_bgd / ds0 * (1-fraction_of_signal)
    # 
    return (sWeights_snl.f, sWeights_bgd.f)
end