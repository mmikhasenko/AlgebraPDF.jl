

@testset "TwoNamedTuples constructors" begin
    @test TwoNamedTuples((a=1.0,)) isa TwoNamedTuples
    @test TwoNamedTuples(; a=1.0) isa TwoNamedTuples
    @test TwoNamedTuples(a=1.0) isa TwoNamedTuples
end