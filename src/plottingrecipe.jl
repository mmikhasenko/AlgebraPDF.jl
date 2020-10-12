
@recipe function f(d::pdf, norm::T where T<:Real=1.0)
    xv = range(d.lims..., length=100)
    return (xv, norm .* d(xv))
end

