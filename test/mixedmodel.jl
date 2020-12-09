
pdf1 = pdf(@. (x;p)->x^2+p.c1, (c1=1.1, ), (-1, 2))
pdf2 = pdf(@. (x;p)->-x^2+p.c2, (c1=1.1, ), (-1, 2))
#
mm′ = MixedModel(SVector(pdf1, pdf2), (f1=0.3,))
mm  = MixedModel([pdf1, pdf2], (f1=0.3,))
#
@test mm(1.1) == mm′(1.1)
@test Ncomp(mm) == 2