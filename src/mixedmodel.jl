
struct MixedModel{T}
    components::SVector{T}
    fractions::NamedTuple
end
#
collectpars(mm::MixedModel) = merge(collectpars.(mm.components)..., mm.fractions)
function fractions(mm::MixedModel; p=mm.fractions)
    updated_fractions = update_namedtuple(mm.fractions,p)
    return vcat(updated_fractions..., (1-sum(updated_fractions)))
end

update_namedtuple(p,from_p) = NamedTuple{Tuple(keys(p))}(getproperty.(Ref(from_p),keys(p)))
# calls
function (mm::MixedModel)(x; p=collectpars(mm))
    _fractions = fractions(mm; p=p)
    v = sum(d(x;p=p)*f for (d,f) in zip(mm.components, _fractions))
    ! prod(iszero, isnan.(v)) &&  @show p#, v
    return v
end
(mm::MixedModel)(x, v) = mm(x; p=v2p(v,mm))

noparsnormf(mm::MixedModel; p=collectpars(mm)) = MixedModel(
        SVector(noparsnormf.(mm.components; p=p)),
        update_namedtuple(mm.fractions,p))