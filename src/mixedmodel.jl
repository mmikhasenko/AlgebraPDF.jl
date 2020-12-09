
struct MixedModel{N} <: AdvancedFunction
    components::SVector{N}
    fractions::NamedTuple
end
# additional constructor
MixedModel(c::Vector, f::NamedTuple) = MixedModel(SVector{length(c)}(c...), f)
# 
Ncomp(mm::MixedModel{N}) where N = N
#
lims(mm::MixedModel) = lims(mm.components[1])
collectpars(mm::MixedModel) = merge(collectpars.(mm.components)..., mm.fractions)
function func(mm::MixedModel,x;p)
    fracs = fractions(mm; p=p)
    v = sum(func(d,x;p=p)*f for (d,f) in zip(mm.components, fracs))
    ! prod(iszero, isnan.(v)) &&  @show p#, v
    return v
end
function integrals(mm::MixedModel, lims; p=collectpars(mm))
    fracs = fractions(mm; p=p)
    ints = integral.(mm.components, Ref(lims); p=p)
    return ints .* fracs
end
integral(mm::MixedModel, lims; p=collectpars(d)) = sum(integrals(mm, lims; p=p))


# inner working
function fractions(mm::MixedModel; p=mm.fractions)
    updated_fractions = updatepars(mm.fractions, p)
    return vcat(updated_fractions..., (1-sum(updated_fractions)))
end
# calls
function (mm::MixedModel)(x; p=collectpars(mm))
    fracs = fractions(mm; p=p)
    v = sum(d(x;p=p)*f for (d,f) in zip(mm.components, fracs))
    ! prod(iszero, isnan.(v)) &&  @show p#, v
    return v
end
(mm::MixedModel)(x, v) = mm(x; p=v2p(v,mm))

noparsnormf(mm::MixedModel; p=collectpars(mm)) = MixedModel(
    SVector(noparsnormf.(mm.components; p=p)),
    updatepars(mm.fractions,p))

# fix parameters
fixpars(mm::MixedModel, pars::NamedTuple) = MixedModel(
    SVector([fixpars(c, selectintersect(collectpars(c), pars)) for c in mm.components]),
    updatepars(mm.fractions, pars))
#
