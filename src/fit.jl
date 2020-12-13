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
cov(fr::FitResults) = invH(fr)
errors(fr::FitResults) = sqrt.(diag(fr.state.invH))
# 
invexacthessian(fr::FitResults, p = minimizer(fr)) = inv(fr.forwarddiff_hessian_callback(p))

function measurements(fr, exacthessian::Bool=false)
    mv = minimizer(fr)
    δv = !(exacthessian) ? errors(fr) : sqrt.(diag(invexacthessian(fr)))
    return [±(m,δ) for (m,δ) in zip(mv, δv)]
end

function fit_llh(data, f; init_pars = error("init_pars!!"), weights = fill(1.0, length(data)))
    llh(p) = -sum((v>0) ? w*log(v) : -1e4 for (w,v) in zip(weights, f(data, p)))
    llh(init_pars) # test
    #
    obj = OnceDifferentiable(llh, init_pars, autodiff = :forwarddiff)
    # 
    m = BFGS()
    options = Optim.Options(show_trace = true)
    lbfgsstate = Optim.initial_state(m, options, obj, init_pars)
    #
    optres = optimize(obj, init_pars, m, options, lbfgsstate)
    #
    hessian_callback(p) = ForwardDiff.hessian(llh, p)
    #
    return FitResults(optres, lbfgsstate, hessian_callback)
end

function fit_llh(data, d::T where T<:AdvancedFunction)
    filtered_data = filter(x->inrange(x, lims(d)), data)
    return fit_llh(filtered_data, d; init_pars=p2v(d))
end
