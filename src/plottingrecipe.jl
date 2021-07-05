
@recipe function f(d::Union{AbstractPDF,MixedModel}, norm::T where T<:Real=1.0, Nbins::Int=100)
    xv = range(lims(d)..., length=Nbins+1)
    return (xv, norm .* d(xv))
end

@recipe f(::Type{T}, d::T) where {T<:AbstractFunctionWithParameters} = x->func(d,x;p=pars(d))

scaletobinneddata(Nd,lims,Nbins) = Nd * (lims[2]-lims[1]) / Nbins
scaletobinneddata(Nd, bins) = Nd * (bins[end]-bins[1]) / (length(bins)-1)


# plotting curve on top of the data
bincenters(bins) = (bins[1:end-1]+bins[2:end]) ./ 2
yerror(y) = sqrt.(y)
bindiffs(x) = x[2:end]- x[1:end-1]


@recipe function f(data::AbstractArray, d::AbstractFunctionWithParameters, Nbins::Integer=60; datalabel="data")
    bins = range(lims(d)..., length=Nbins)
    @series begin
        centers = bincenters(bins)
        bincontent = [sum(x->l<x<r, data) for (l,r) in zip(bins[1:end-1], bins[2:end])]
        # 
        x := centers
        y := bincontent
        xerror := bindiffs(centers) / 2
        yerror := yerror(bincontent)
        seriestype := :scatter
        seriescolor --> :black
        markersize --> 3
        label --> datalabel
        ()
    end
    @series begin
        (d, scaletobinneddata(length(data), bins))
    end
end

# two dimensional
const ArrayOrRange = Union{AbstractArray,AbstractRange}
# 
@recipe f(
    xv::ArrayOrRange, yv::ArrayOrRange,
    d::AbstractFunctionWithParameters) = (xv,yv,(x,y)->func(d,(x,y)))
