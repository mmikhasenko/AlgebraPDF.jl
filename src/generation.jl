struct binned1dDensity{T}
    grid::Vector{Float64}
    cumarr::Vector{T}
    density::Function
end

function binned1dDensity(grid, weights; density::Function = x->(println("warning: no callback provided!");1))
    _weights = weights ./ sum(weights)
    cumarr = [0;cumsum(vcat(_weights...), dims=1)]
    return binned1dDensity(grid, cumarr, density)
end

function getbinned1dDensity(g, lims, Nbins)
    grid = collect(range(lims[1], lims[2], length=Nbins))
    #
    weights = [g(v) for v in (grid[2:end] .+ grid[1:end-1]) ./ 2]
    return binned1dDensity(grid, weights; density=g)
end

import Base:rand
function rand(bD::binned1dDensity; checknonzero::Bool=false)
    # get
    binind = findfirst(bD.cumarr .> rand())-1
    σl, σr = bD.grid[binind], bD.grid[binind+1]
    σ = σl+rand()*(σr-σl)
    checknonzero && bD.density(σ) == 0.0 && return rand(bD,checknonzero)
    return σ
end

function generate(Nev::Int, d::T where T <: AdvancedFunction; p=collectpars(d), Nbins=100)
    grid = collect(range(lims(d)..., length=Nbins))
    centers = (grid[2:end] .+ grid[1:end-1]) ./ 2
    #
    weights = d(centers)
    bD = binned1dDensity(grid, weights; density=x->d(x))
    return [rand(bD) for _ in 1:Nev]
end
