using AlgebraPDF # https://github.com/mmikhasenko/AlgebraPDF.j
using BenchmarkTools

const αv = 0.37
const βv = -0.002
const σev = [0.58, 1.19, 1.35, 1.58] .* 1e-3;
const mth = 0.96;
# 
const Mv = [0.99, 1.02, 1.06, 1.11]
const Γv = [0.005, 0.001, 0.004, 0.010]
# 
const fitlims = (0, 0.22);

pq(s) = s< mth^2 ? 0.0 : sqrt(s-mth^2)/sqrt(s)
BW(s,m,Γ) = 1/(m^2-s-1im*m*Γ)
BWmth(e,ΔM,Γ) = (e+mth)*Γ*BW((e+mth)^2, ΔM+mth, Γ)
pqmth(e) = pq((e+mth)^2)
bgd(e,α,β) = e^α*exp(-e*β)

# construction of PDFs
function construct_pdf_directly_integral_conv()
    function four_peaks_and_background(e; p)
        σev = [0.58, 1.19, 1.35, 1.58] .* 1e-3
             (conv_with_gauss.(e, y->abs2(BWmth(y, p.m1, p.Γ1)), σev[1]) .+
        p.f2.*conv_with_gauss.(e, y->abs2(BWmth(y, p.m2, p.Γ2)), σev[2]) .+
        p.f3.*conv_with_gauss.(e, y->abs2(BWmth(y, p.m3, p.Γ3)), σev[3]) .+
        p.f4.*conv_with_gauss.(e, y->abs2(BWmth(y, p.m4, p.Γ4)), σev[4]) .+
        p.fb.*bgd.(e, p.α, p.β)).*pqmth.(e)
    end
    #
    pdf_direct = 
        Normalized(four_peaks_and_background; lims=fitlims, p=(
            m1=Mv[1]-mth, Γ1=Γv[1],
            m2=Mv[2]-mth, Γ2=Γv[2],
            m3=Mv[3]-mth, Γ3=Γv[3],
            m4=Mv[4]-mth, Γ4=Γv[4],
            f2=5.0, f3=0.2, f4=0.1, fb=0.1,
            α = αv, β = βv
        ))
    return pdf_direct
end

function construct_pdf_modular_integral_conv()
    pdf1 = Normalized(@. (e;p)->abs2(BWmth(e, p.m1, p.Γ1));  p = (m1=Mv[1]-mth, Γ1=Γv[1]), lims = fitlims)
    pdf2 = Normalized(@. (e;p)->abs2(BWmth(e, p.m2, p.Γ2));  p = (m2=Mv[2]-mth, Γ2=Γv[2]), lims = fitlims)
    pdf3 = Normalized(@. (e;p)->abs2(BWmth(e, p.m3, p.Γ3));  p = (m3=Mv[3]-mth, Γ3=Γv[3]), lims = fitlims)
    pdf4 = Normalized(@. (e;p)->abs2(BWmth(e, p.m4, p.Γ4));  p = (m4=Mv[4]-mth, Γ4=Γv[4]), lims = fitlims)
    pdfb = Normalized(@. (e;p)->bgd(e, αv, βv);  p = NamedTuple(), lims = fitlims)
    pdfS = [
        conv_with_gauss(pdf1, σev[1]),
        conv_with_gauss(pdf2, σev[2]),
        conv_with_gauss(pdf3, σev[3]),
        conv_with_gauss(pdf4, σev[4]),
        pdfb]
    #
    phsp = Normalized(@. (e;p)->pqmth(e);  p = NamedTuple(), lims = fitlims)
    pdfS = [d * phsp for d in pdfS]
    #
    pdfS[2] *= (f2=5.0,)
    pdfS[3] *= (f3=0.2,)
    pdfS[4] *= (f4=0.1,)
    pdfS[5] *= (fb=0.1,)
    # 
    pdf_modular = sum(pdfS);
    return pdf_modular
end



function construct_pdf_light_modular_integral_conv()
    pdf1 = Normalized(@. (e;p)->abs2(BWmth(e, p.m1, p.Γ1)) * pqmth(e) * 1.0 ;  p = (m1=Mv[1]-mth, Γ1=Γv[1]          ), lims = fitlims)
    pdf2 = Normalized(@. (e;p)->abs2(BWmth(e, p.m2, p.Γ2)) * pqmth(e) * p.f2;  p = (m2=Mv[2]-mth, Γ2=Γv[2], f2 = 5.0), lims = fitlims)
    pdf3 = Normalized(@. (e;p)->abs2(BWmth(e, p.m3, p.Γ3)) * pqmth(e) * p.f3;  p = (m3=Mv[3]-mth, Γ3=Γv[3], f3 = 0.2), lims = fitlims)
    pdf4 = Normalized(@. (e;p)->abs2(BWmth(e, p.m4, p.Γ4)) * pqmth(e) * p.f4;  p = (m4=Mv[4]-mth, Γ4=Γv[4], f4 = 0.1), lims = fitlims)
    pdfb = Normalized(@. (e;p)->bgd(e, αv, βv) * pqmth(e) * p.f4;  p = (fb=0.1,), lims = fitlims)
    pdfS = [
        conv_with_gauss(pdf1, σev[1]),
        conv_with_gauss(pdf2, σev[2]),
        conv_with_gauss(pdf3, σev[3]),
        conv_with_gauss(pdf4, σev[4]),
        pdfb]
    #
    return pdfS
end

# benchmark
const xv = range(fitlims..., length=10)
# 
pdf_direct = construct_pdf_directly_integral_conv()
@btime pdf_direct(0.1) # 19 ms
@btime pdf_direct(xv)  # 20 ms
#
pdf_modular = construct_pdf_modular_integral_conv()
@btime pdf_modular(0.1) # 182 ms
@btime pdf_modular(xv) # 190 ms 
#
pdfs_light = construct_pdf_light_modular_integral_conv()
pdf_light = sum(pdfs_light)
@btime pdf_light(0.1) # 225 ms
@btime pdf_modular(xv) # 197 ms 

using Plots
using HEPRecipes # https://github.com/mmikhasenko/HEPRecipes.jl
let 
    Nev = 1000
    data = generate(Nev, pdf_direct; Nbins=500)
    # 
    bins = range(fitlims..., length=100)
    N = Nev / (length(bins)-1) * (fitlims[2]-fitlims[1])
    #
    fine = range(fitlims..., length=300)
    plot(fine, N .* pdf_direct(fine), l=(2, :red), lab="model")
    #
    wh = weightedHistogram(data, bins=bins)
    plot!(wh, m=(3,:black), l=(2, :red), lab="data")
end
