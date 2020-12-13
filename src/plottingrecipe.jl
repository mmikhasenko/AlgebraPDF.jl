
@recipe function f(d::T where T<:AdvancedFunction, norm::T where T<:Real=1.0, Nbins::Int=100)
    xv = range(lims(d)..., length=Nbins+1)
    return (xv, norm .* d(xv))
end

scaletobinneddata(Nd,lims,Nbins) = Nd * (lims[2]-lims[1]) / Nbins
scaletobinneddata(Nd, bins) = Nd * (bins[end]-bins[1]) / (length(bins)-1)
