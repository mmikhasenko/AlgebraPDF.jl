using Plots
using AlgebraPDF
using Test

let
    p = aGauss((μ=1.1, σ = 2.2), (0,5))
    data = generate(1000, p)
    plot(data, p, lab="m1")
end