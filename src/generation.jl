struct binned1dDensity{T}
    grid::Vector{Float64}
    cumarr::Vector{T}
    density::Function
end

function binned1dDensity(
    grid,
    weights;
    density::Function = x -> (println("warning: no callback provided!"); 1),
)
    _weights = weights ./ sum(weights)
    cumarr = [0; cumsum(vcat(_weights...), dims = 1)]
    return binned1dDensity(grid, cumarr, density)
end

function getbinned1dDensity(g, lims, Nbins)
    grid = collect(range(lims[1], lims[2], length = Nbins))
    #
    weights = [g(v) for v in (grid[2:end] .+ grid[1:end-1]) ./ 2]
    return binned1dDensity(grid, weights; density = g)
end

import Base: rand
function rand(bD::binned1dDensity; removenegative::Bool = false)
    # get
    binind = findfirst(bD.cumarr .> rand()) - 1
    σl, σr = bD.grid[binind], bD.grid[binind+1]
    σ = σl + rand() * (σr - σl)
    removenegative && bD.density(σ) == 0.0 && return rand(bD, removenegative)
    return σ
end

function generate(
    d::T where {T<:AbstractFunctionWithParameters},
    N::Int;
    p = freepars(d),
    Nbins = 100,
    removenegative::Bool = false,
)
    grid = collect(range(lims(d)..., length = Nbins))
    centers = (grid[2:end] .+ grid[1:end-1]) ./ 2
    #
    weights = d(centers; p)
    bD = binned1dDensity(grid, weights; density = x -> d(x; p))
    return [rand(bD; removenegative) for _ = 1:N]
end
rand(d::T where {T<:AbstractFunctionWithParameters}, N::Int = 1) = generate(d, N)
