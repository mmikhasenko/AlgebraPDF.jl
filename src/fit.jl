#
llh(data, f; weights = fill(1.0, length(data))) = -sum((v>0) ? w*log(v) : -1e4 for (w,v) in zip(weights, f(data)))

struct FitResults{T<:Optim.OptimizationResults, V<:Optim.AbstractOptimizerState} 
    optres::T
    state::V
    forwarddiff_hessian_callback::Function
end
minimizer(fr::FitResults) = minimizer(fr.optres)
minimum(fr::FitResults) = minimum(fr.optres)
# 
invH(fr::FitResults) = fr.state.invH
covmat(fr::FitResults) = invH(fr)
errors(fr::FitResults) = sqrt.(diag(fr.state.invH))
# 
invexacthessian(fr::FitResults, p = minimizer(fr)) = inv(fr.forwarddiff_hessian_callback(p))

function measurements(fr, exacthessian::Bool=false)
    mv = minimizer(fr)
    δv = !(exacthessian) ? errors(fr) : sqrt.(diag(invexacthessian(fr)))
    return [±(m,δ) for (m,δ) in zip(mv, δv)]
end

#      _|_|  _|    _|      
#    _|          _|_|_|_|  
#  _|_|_|_|  _|    _|      
#    _|      _|    _|      
#    _|      _|      _|_|  


function optimize_get_state_hessian(func, init_pars, m = BFGS(); show_trace::Bool = true)
    obj = OnceDifferentiable(func, init_pars, autodiff = :forwarddiff)
    # 
    options = Optim.Options(show_trace = show_trace)
    lbfgsstate = Optim.initial_state(m, options, obj, init_pars)
    #
    optres = optimize(obj, init_pars, m, options, lbfgsstate)
    #
    hessian_callback(p) = ForwardDiff.hessian(func, p)
    #
    return FitResults(optres, lbfgsstate, hessian_callback)
end

function fit_llh(data, f; init_pars = error("init_pars!!"), weights = fill(1.0, length(data)))
    llh(p) = -sum((v>0) ? w*log(v) : -1e4 for (w,v) in zip(weights, f(data, p)))
    return optimize_get_state_hessian(llh, init_pars, BFGS())
end

function fit_llh_with_constraints(data, f, fc; init_pars = error("init_pars!!"))
    llh(p) = -sum((v>0) ? log(v) : -1e4 for v in f(data, p))
    llh_constr(p) = llh(p) + fc(p)
    return optimize_get_state_hessian(llh_constr, init_pars, BFGS())
end



# AdvancedFunction
function fit_llh(data, d::T where T<:AdvancedFunction)
    filtered_data = filter(x->inrange(x, lims(d)), data)
    return fit_llh(filtered_data, d; init_pars=p2v(d))
end

chi2(specification; p) = sum(((getproperty(p, k)-v) / e)^2
    for (k,(v,e)) in zip(keys(specification),specification))

function fit_llh_with_constraints(data, d::T where T<:AdvancedFunction, specification)
    filtered_data = filter(x->inrange(x, lims(d)), data)
    fc(v) = chi2(specification; p=v2p(v, d))
    return fit_llh_with_constraints(filtered_data, d, fc; init_pars=p2v(d))
end
