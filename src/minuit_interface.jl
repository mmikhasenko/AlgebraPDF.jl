# IMinuit initialization
const iminuit = PyNULL()

function load_python_deps!()
    copy!(iminuit, pyimport_conda("iminuit", "iminuit", "conda-forge"))
    return nothing
end

function __init__()
    load_python_deps!()
    # try
    # catch ee
    #     if PyCall.conda
    #         Conda.add("iminuit")
    #         load_python_deps!()
    #     else
    #         typeof(ee) <: PyCall.PyError || rethrow(ee)
    #         @warn("""
    #              Python Dependencies not installed!
    #              """)
    #     end
    # end
    # return nothing
end

#######################################################################

abstract type AbstractMinuit end
struct MigradAndHesse <: AbstractMinuit
    errordef::Float64
end
MigradAndHesse() = MigradAndHesse(0.5)  # 1/2 for likelihood

#######################################################################

@with_kw struct MinuitAndHesseResults{T,N,R}
    fmin::T
    params::N
    hesse::R
end
#
const keys_FMin = (:fval,
    :edm,
    :tolerance,
    :nfcn,
    :nfcn_total,
    :up,
    :is_valid,
    :has_valid_parameters,
    :has_accurate_covar,
    :has_posdef_covar,
    :has_made_posdef_covar,
    :hesse_failed,
    :has_covariance,
    :is_above_max_edm,
    :has_reached_call_limit,
    :has_parameters_at_limit,
    :ngrad,
    :ngrad_total)

const keys_Param = (:number,
    :name,
    :value,
    :error,
    :is_const,
    :is_fixed,
    :has_limits,
    :has_lower_limit,
    :has_upper_limit,
    :lower_limit,
    :upper_limit)
#

const str_class_MigradResult = "PyObject <class 'iminuit.util.MigradResult'>"
const str_class_Params = "PyObject <class 'iminuit.util.Params'>"

#
"""
    nt(obj::PyObject)

Function to parse the Minuit output to julia named tuples
"""
function nt(obj::PyObject)
    if string(obj.__class__) == str_class_MigradResult
        fmin = NamedTuple{keys_FMin}(obj.fmin)
        params = NamedTuple{keys_Param}.(obj.params)
        return (fmin = fmin, params = params)
    end
    if string(obj.__class__) == str_class_Params
        return NamedTuple{keys_Param}.(obj)
    end
    error("unknown PyObject")
end

#######################################################################


function minimize(fcn, init_pars, optimizer::MigradAndHesse; kws...)
    # 
    o = iminuit.Minuit.from_array_func(fcn, collect(init_pars);
        pedantic = false, errordef=optimizer.errordef)
    println("Minuit object is created: fit with $(length(init_pars)) parameters")
    # 
    migrad_result = pycall(o.migrad, PyObject)
    println("MIGRAD is finishied")
    # 
    hesse_result = pycall(o.hesse, PyObject)
    println("HESSE is finishied")
    # 
    # 
    MinuitAndHesseResults(; hesse = nt(hesse_result), nt(migrad_result)...)
end

parameternames(mr::MinuitAndHesseResults) = Tuple(Symbol.(getproperty.(mr.params, :name)))
minimizer(mr::MinuitAndHesseResults) = getproperty.(mr.params, :value)
errors(mr::MinuitAndHesseResults) = getproperty.(mr.hesse, :error)
measurements(mr::MinuitAndHesseResults) = map(x->±(x...), zip(minimizer(mr), errors(mr)))
minimum(mr::MinuitAndHesseResults) = mr.fmin.fval

