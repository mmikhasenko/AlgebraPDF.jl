# #######################################################################
# Minuit MWE
# o = AlgebraPDF.iminuit.Minuit(x->x[1]^2+x[2]^2, [1.0, 1.0])
# o.errordef = 1/2
# migrad_result = AlgebraPDF.pycall(o.migrad, AlgebraPDF.PyObject)
# hesse_result = AlgebraPDF.pycall(o.hesse, AlgebraPDF.PyObject)
# 
# migrad_result.fmin
# migrad_result.params
# migrad_result.covariance

# IMinuit initialization
const iminuit = PyNULL()

function __init__()
    copy!(iminuit, pyimport_conda("iminuit", "iminuit=2.7.0", "conda-forge"))
    iminuit.__version__ < "2.7.0" && println("""
    Version of iminuit >= v2.7.0 is required, try to install it manually using `PyCall.Conda`:

    julia> using PyCall
    julia> PyCall.Conda.add("iminuit=2.7.0")

    check the version with

    julia> iminuit = pyimport_conda("iminuit", "iminuit=2.7.0", "conda-forge")
    julia> iminuit.__version__

    or 

    julia> iminuit.__version__
    """)
end

isiminuitimported() = (AlgebraPDF.iminuit != AlgebraPDF.PyNULL())

#######################################################################

abstract type AbstractMinuit end
struct MigradAndHesse <: AbstractMinuit
    errordef::Float64
end
MigradAndHesse() = MigradAndHesse(0.5)  # 1/2 for likelihood

# #######################################################################
@with_kw struct MinuitAndHesseResults{T,N,R}
    fmin::T
    params::N
    hesse::R
end
#
const keys_FMin = (:algorithm,
    :edm,
    :edm_goal,
    :errordef,
    :fval,
    :has_accurate_covar,
    :has_covariance,
    :has_made_posdef_covar,
    :has_parameters_at_limit,
    :has_posdef_covar,
    :has_reached_call_limit,
    :has_valid_parameters,
    :hesse_failed,
    :is_above_max_edm,
    :is_valid,
    :nfcn,
    :ngrad,
    :reduced_chi2)
# 

const keys_Param = (
    :number,
    :name,
    :value,
    :error,
    :merror,
    :is_const,
    :is_fixed,
    :lower_limit,
    :upper_limit)
#

const str_class_Minuit = "PyObject <class 'iminuit.minuit.Minuit'>"
const str_class_Params = "PyObject <class 'iminuit.util.Param'>"

obj2nt(obj, keys) = NamedTuple{keys}(getproperty(obj, s) for s in keys)
function obj2nt(obj::PyObject)
    if string(obj.__class__) == str_class_Minuit
        fmin = obj2nt(obj.fmin, AlgebraPDF.keys_FMin)
        params = obj2nt.(obj.params)
        return (fmin = fmin, params = params)
    end
    if string(obj.__class__) == str_class_Params
        return obj2nt(obj, keys_Param)
    end
    error("unknown PyObject")
end

#######################################################################

function minimize(fcn, init_pars, optimizer::MigradAndHesse; kws...)
    # 
    o = iminuit.Minuit(fcn, collect(init_pars))
    o.errordef = optimizer.errordef
    println("Minuit object is created: fit with $(length(init_pars)) parameters")
    # 
    migrad_result = pycall(o.migrad, PyObject)
    println("MIGRAD is finishied")
    # 
    hesse_result = pycall(o.hesse, PyObject)
    println("HESSE is finishied")
    # 
    MinuitAndHesseResults(; hesse = hesse_result.covariance, obj2nt(migrad_result)...)
end

parameternames(mr::MinuitAndHesseResults) = Tuple(Symbol.(getproperty.(mr.params, :name)))
minimizer(mr::MinuitAndHesseResults) = getproperty.(mr.params, :value)
errors(mr::MinuitAndHesseResults) = getproperty.(mr.params, :error)
measurements(mr::MinuitAndHesseResults) = map(x->±(x...), zip(minimizer(mr), errors(mr)))
minimum(mr::MinuitAndHesseResults) = mr.fmin.fval

