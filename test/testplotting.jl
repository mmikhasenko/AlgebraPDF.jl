using AlgebraPDF
using Test

@testset "Plotting utils" begin
    @test scaletobinneddata(10, (0,1),10) ≈ 1.0
    @test scaletobinneddata(10, range(0,1,length=11)) ≈ 1.0
end
