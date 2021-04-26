
struct MixedModel{N1,N2,T} <: AbstractFunctionWithParameters
    components::SVector{N1}
    fractions::T
    keys::SVector{N2,Symbol}
end
# additional constructor
MixedModel(c::Vector, f::NamedTuple) = MixedModel(SVector{length(c)}(c...), Parameters(f), SVector{length(c)-1}(keys(f)...))

function nt_fractions(N)
    f = 1/N
    fs = [Symbol("f",i) for i in 0:N-2]
    NamedTuple{Tuple(fs)}(Iterators.repeat([f], N-1))
end
function MixedModel(c::Vector)
    fractions = nt_fractions(length(c))
    println("
    List of fraction parameters is created automatically
       $(fractions)
       Make sure that it does not overlar with parameters of the component")
    MixedModel(c, nt_fractions(length(c)))
end
# 
Ncomp(mm::MixedModel{N}) where N = N
fractions(mm::MixedModel) = mm.fractions
#
lims(mm::MixedModel) = lims(mm.components[1])

pars(mm::MixedModel) = sum(pars.(mm.components)) + fractions(mm)
function func(mm::MixedModel,x;p)
    fracs = fractionvalues(mm; p=p)
    v = sum(func(d,x;p=p)*f for (d,f) in zip(mm.components, fracs))
    ! prod(iszero, isnan.(v)) &&  @show p#, v
    return v
end

# inner working
orderedfractionvalues(fs::Parameters, keys) = [getproperty(fs, k) for k in keys]
# 
function fractionvalues(mm::MixedModel; p=∅)
    updated_fractions = updatepars(mm.fractions, p)
    ordered_fractions= orderedfractionvalues(updated_fractions, mm.keys)
    return vcat(ordered_fractions..., (1-sum(ordered_fractions)))
end
# calls
function (mm::MixedModel)(x; p=freepars(mm))
    fracs = fractionvalues(mm; p=p)
    v = sum(d(x;p=p)*f for (d,f) in zip(mm.components, fracs))
    ! prod(iszero, isnan.(v)) &&  @show p#, v
    return v
end
(mm::MixedModel)(x, v) = mm(x; p=v2p(v,mm))

noparsnormf(mm::MixedModel; p=freepars(mm)) = MixedModel(
    SVector(noparsnormf.(mm.components; p=p)),
    updatepars(mm.fractions,p))

# fix parameters
copy(mm::MixedModel, p) = MixedModel(
    SVector([copy(c, p) for c in mm.components]),
    copy(mm.fractions, p), mm.keys)

# 
"""
    integrals(mm::MixedModel, lims; p=freepars(mm))

Computes integrals of the components in a given range.
"""
function integrals(mm::MixedModel, lims; p=freepars(mm))
    fracs = fractionvalues(mm; p=p)
    ints = integral.(mm.components, Ref(lims); p=p)
    return ints .* fracs
end
integral(mm::MixedModel, lims; p=freepars(mm)) = sum(integrals(mm, lims; p=p))
