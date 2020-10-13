struct xProductPDF{N}
    dims::SVector{N,pdf}
    keys::SVector{N,Symbol}
end

xProductPDF(; kwargs...) = xProductPDF(SVector(values(kwargs)...), SVector(keys(kwargs)...))

function generate(N, X::xProductPDF)
    data_vectors = generate.(N, X.dims)
    ntv = NamedTuple{Tuple(X.keys)}.(zip(data_vectors...))
    return ntv
end