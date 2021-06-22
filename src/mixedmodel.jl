
const SequenceType = Union{Vector, Tuple, NamedTuple, SVector}

struct MixedModel{N,T} <: AbstractFunctionWithParameters
    components::SVector{N}
    fractions::T
end
# additional constructor
MixedModel(c::SequenceType, f) = MixedModel(SVector{length(c)}(c...), f)

fractions(mm::MixedModel) = mm.fractions


# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # #
# # to be implemented for F <: AbstractFunctionWithParameters
# pars(d::F, isfree::Bool)
# func(d::F, x; p)
# updatevalueorflag( d::F, s::Symbol, isfree::Bool, v=getproperty(pars(d),s))
# function (mm::F)(x; p=freepars(mm))
# # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # # 


pars(mm::MixedModel, isfree::Bool) = sum(pars.(mm.components, isfree)) + pars(fractions(mm), isfree)

function func(mm::MixedModel,x;p)
    fracs = fractionvalues(mm; p=p)
    v = sum(func(d,x;p=p)*f for (d,f) in zip(mm.components, fracs))
    ! prod(iszero, isnan.(v)) &&  @show p#, v
    return v
end

function updatevalueorflag(d::MixedModel{N,T}, s::Symbol, isfree::Bool, v=getproperty(pars(d),s)) where {N,T}
    newcomponents = []
    for c in d.components
        push!(newcomponents,
            ispar(c,s) ? updatevalueorflag(c, s, isfree, v) : c)
    end
    #
    fs = fractions(d)
    return MixedModel(
        SVector{N}(newcomponents),
        ispar(fs,s) ? updatevalueorflag(fs, s, isfree, v) : fs)
end

# calls
function (mm::MixedModel)(x; p=freepars(mm))
    fracs = fractionvalues(mm; p=p)
    v = sum(d(x; p=p)*f for (d,f) in zip(mm.components, fracs))
    ! prod(iszero, isnan.(v)) &&  @show p#, v
    return v
end
(mm::MixedModel)(x, v) = mm(x; p=v2p(v,mm))




# extra functions

function defaultfractionsnt(N, fraction_letter)
    names = [Symbol(fraction_letter,i) for i in 1:N-1]
    values = Iterators.repeat([1/N], N-1)
    return NamedTuple{Tuple(names)}(values)
end

"""
    MixedModel(c::Vector; fraction_letter="f")

Constructor that creates the fraction named tuple automatically, (f1=1/N, f2=1/N, ..., fN-1=1/N)
"""
function MixedModel(c::Vector; fraction_symbol="f")
    fractions = defaultfractionsnt(length(c), fraction_symbol)
    println("
    List of fraction parameters is created automatically
       $(fractions)
       Make sure that it does not overlap with parameters of the component")
    return MixedModel(c, fractions)
end

#
Ncomp(mm::MixedModel{N}) where N = N
#
lims(mm::MixedModel) = lims(mm.components[1])


# inner working
fractionvalues(mm::MixedModel; p=∅) = vcat(mm.fractions..., (1-sum(mm.fractions)))


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


noparsnormf(mm::MixedModel; p=freepars(mm)) = MixedModel(
    SVector(noparsnormf.(mm.components; p=p)),
    updatepars(mm.fractions,p))
