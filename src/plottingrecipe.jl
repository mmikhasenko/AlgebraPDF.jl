
@recipe function f(d::Union{AbstractPDF,MixedModel,SumOfPDF}, norm::T where T<:Real=1.0, bins::Int=100)
    xv = range(lims(d)..., length=bins+1)
    return (xv, norm .* d(xv))
end

@recipe f(::Type{T}, d::T) where {T<:AbstractFunctionWithParameters} = x->func(d,x;p=pars(d))



scaletobinneddata(Nd,lims,Nbins) = Nd * (lims[2]-lims[1]) / Nbins
scaletobinneddata(Nd, bins) = Nd * (bins[end]-bins[1]) / (length(bins)-1)


# plotting curve on top of the data
bincenters(bins) = (bins[1:end-1]+bins[2:end]) ./ 2
yerror(y) = sqrt.(y)
bindiffs(x) = x[2:end]- x[1:end-1]


@recipe function f(data::AbstractArray, d::Union{AbstractPDF,MixedModel}; bins=60, datalabel="data")
    binning = range(lims(d)..., length=bins+1)
    @series begin
        (d, scaletobinneddata(length(data), binning))
    end
    @series begin
        centers = bincenters(binning)
        bincontent = [sum(x->l<x<r, data) for (l,r) in zip(binning[1:end-1], binning[2:end])]
        # 
        x := centers
        y := bincontent
        xerror := bindiffs(centers) / 2
        yerror := yerror(bincontent)
        seriestype := :scatter
        seriescolor --> :black
        markersize --> 3
        label := datalabel
        ()
    end
end

# 
@recipe f(
    xv::ArrayOrRange, yv::ArrayOrRange,
    d::AbstractFunctionWithParameters) = (xv,yv,(x,y)->func(d,(x,y)))
