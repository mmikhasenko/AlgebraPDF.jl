using AlgebraPDF
using Test

@makefuntype Amazing2D(x;p) = (x[1]-p.x0)^2+(x[2]-p.y0)^2-p.R0^2
a = Amazing2D((x0=1.1, y0=2.1, R0=0.0))

@testset "Two-dimensional function" begin
    @test a((1.1,2.1)) == 0
    @test a((1.1+1,2.1+1)) ≈ 2
    data = collect(zip(rand(10), rand(10)))
    @test length(a(data)) == 10
end


@makefuntype Amazing3D(x;p) = (x[1]-p.x0)^2+(x[2]-p.y0)^2+(x[3]-p.z0)^2-p.R0^2
a = Amazing3D((x0=1.1, y0=2.1, z0=3.1, R0=0.0))

@testset "Two-dimensional function" begin
    @test a((1.1,2.1,3.1)) == 0
    @test a((1.1+1,2.1+1,3.1+1)) ≈ 3
    data = collect(zip(rand(10), rand(10), rand(10)))
    @test length(a(data)) == 10
end


