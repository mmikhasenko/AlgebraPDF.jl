struct convGauss{T,S} <: AbstractPDF{1}
    source::T
    σ::S
end

pars(d::convGauss, isfree::Bool) = merge(pars(d.source, isfree), pars(d.σ, isfree))
lims(d::convGauss) = lims(d.source)
updatevalueorflag(d::convGauss, s::Symbol, isfree::Bool, v = getproperty(pars(d), s)) =
    convGauss(
        ispar(d.source, s) ? updatevalueorflag(d.source, s, isfree, v) : d.source,
        ispar(d.σ, s) ? updatevalueorflag(d.σ, s, isfree, v) : d.σ,
    )

function func(d::convGauss, x::NumberOrTuple; p = pars(d))
    σ = func(d.σ, x; p)
    g(z) = AlgebraPDF.standardgauss(z, σ)
    f(z) = func(d.source, z; p)
    return quadgk(y -> f(x - y) * g(y), -5 * σ, +5 * σ)[1]
end
