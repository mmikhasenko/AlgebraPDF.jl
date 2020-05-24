using Test
using AlgebraPDF

BW(s, m, Γ) = 1 / (m^2 - s^2 - 1im*m*Γ)
#
pdf1 = pdf(@. (e;p)->abs2(BW(e^2, p.m1, p.Γ1));
    p0 = (m1=0.25, Γ1=2e-3), lims = (0, 0.15))
#
@test pdf1(rand()) != 0.0
@test length(pdf1(rand(10))) == 10
#
#
pdf2 = pdf(@. (e;p)->abs2(BW(e^2, p.m2, p.Γ2));
    p0 = (m2=0.1, Γ2=14e-3), lims = (0, 0.15))
#
let v = rand(2), x = 1.1
    @test pdf2(1.1, v) == pdf2(1.1; p=v2p(v, pdf2))
end
#
#
pdf2 *= (f2=3.0,)
@test length(pdf2.p0) == 3
#
#
pdf_sum = pdf1 + pdf2
@test npars(pdf_sum) == 5
#
