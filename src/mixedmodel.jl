
struct MixedModel{T}
    components::SVector{T}
    fractions::NamedTuple
end
#
collectpars(mm::MixedModel) = merge(collectpars.(mm.components)..., mm.fractions)

update_namedtuple(p,from_p) = NamedTuple{Tuple(keys(p))}(getproperty.(Ref(from_p),keys(p)))
# calls
function (mm::MixedModel)(x; p=collectpars(mm))
    # @show p
    updated_fractions = update_namedtuple(mm.fractions,p)
    fractions = vcat(updated_fractions..., (1-sum(updated_fractions)))
    v = sum(d(x;p=p)*f for (d,f) in zip(mm.components, fractions))
    ! prod(iszero, isnan.(v)) &&  @show p#, v
    return v
end
(mm::MixedModel)(x, v) = mm(x; p=v2p(v,mm))
