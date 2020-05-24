struct binned1dDensity
    grid::Array{Real,1}
    cumarr::Array{Real}
    density
end

function getbinned1dDensity(g, lims, Nbins)
    grid = [σi for σi=range(lims[1], lims[2], length=Nbins)]
    #
    weights = [g(v) for v in (grid[2:end] .+ grid[1:end-1]) ./ 2]
    weights ./= sum(weights)
    cumarr = [0;cumsum(vcat(weights...), dims=1)]
    return binned1dDensity(grid,cumarr,g)
end

import Base:rand
function rand(bD::binned1dDensity)
    # get
    binind = findfirst(bD.cumarr .> rand())-1
    σl, σr = bD.grid[binind], bD.grid[binind+1]
    σ = σl+rand()*(σr-σl)
    bD.density(σ) == 0.0 && return rand(bD)
    return σ
end

function generate(Nev::Int, d::pdf; p=d.p0, Nbins=100)
    bD = getbinned1dDensity(x->d.f(x; p=p), d.lims, Nbins)
    return [rand(bD) for _ in 1:Nev]
end
