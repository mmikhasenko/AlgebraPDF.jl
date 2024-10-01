using AlgebraPDF
using Test

@testset "Generate in range" begin
    d = Normalized(FGauss((a=0,b=1.1)), (2,5))
    data = generate(d, 100)
    sum(x->inrange(x,lims(d)), data) == length(data)
end