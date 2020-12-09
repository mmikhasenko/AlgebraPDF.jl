
@recipe function f(d::T where T<:AdvancedFunction, norm::T where T<:Real=1.0)
    xv = range(lims(d)..., length=100)
    return (xv, norm .* d(xv))
end

