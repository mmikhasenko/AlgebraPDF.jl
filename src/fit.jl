
function fit(objective::Union{NegativeLogLikelihood, Extended, ChiSq},
        optimizer = MigradAndHesse(), args...; kws...)
    d = model(objective)
    fit_result = minimize(x->objective(0,x), p2v(d), optimizer, args...; kws...)
    return fit_summary(fit_result, d)
end

fit(model::Normalized, data, optimizer = MigradAndHesse(), args...; kws...) =
    fit(NegativeLogLikelihood(model, data), optimizer, args...; kws...)
#

# the optimizer interface <: Optimizer should implement
#  - minimize(objective, optimizer::Optimizer)
#  - minimizer(minimizationResult)
#  - measurement(minimizationResult)
#  - minimum(minimizationResult)


function fit_summary(fit_result, model)
    nll = minimum(fit_result)
    parameters = v2p(minimizer(fit_result), model)
    meas = v2p(measurements(fit_result), model)
    best_model = updatepars(model, parameters)
    #
    return (; parameters, best_model, nll, fit_result, measurements=meas)
end


# function fit(data::AbstractArray, model::AbstractFunctionWithParameters, args...; kws...)
#     fit_results = fit_llh(data, model, args...; init_pars = p2v(model), kws...)
#     return fit_summary(fit_results, model) + (fit_results = fit_results, )
# end








                                                                                                         
#                                            _|                          _|              _|                
#    _|_|_|    _|_|    _|_|_|      _|_|_|  _|_|_|_|  _|  _|_|    _|_|_|      _|_|_|    _|_|_|_|    _|_|_|  
#  _|        _|    _|  _|    _|  _|_|        _|      _|_|      _|    _|  _|  _|    _|    _|      _|_|      
#  _|        _|    _|  _|    _|      _|_|    _|      _|        _|    _|  _|  _|    _|    _|          _|_|  
#    _|_|_|    _|_|    _|    _|  _|_|_|        _|_|  _|          _|_|_|  _|  _|    _|      _|_|  _|_|_|    
                                                                                                         





# function fit_llh(data, d::AbstractFunctionWithParameters, constraints::NamedTuple;
#     init_pars = p2v(d), kws...)
#     # 
#     filtered_data = filter(x->inrange(x, AlgebraPDF.lims(d)), data)
#     # 
#     llh(p) = -sum((v>0) ? log(v) : -1e4 for v in d(filtered_data, p))
#     fc(v) = AlgebraPDF.chi2(constraints; p=v2p(v, d))
#     llh_constr(p) = llh(p) + fc(p)
#     return minimize(llh_constr, init_pars; kws...)
# end


# chi2(specification; p) = sum(((getproperty(p, k)-v) / e)^2
#     for (k,(v,e)) in zip(keys(specification),specification))

# function fit_llh_with_constraints(data, d::AbstractPDF, specification)
#     filtered_data = filter(x->inrange(x, lims(d)), data)
#     fc(v) = chi2(specification; p=v2p(v, d))
#     return fit_llh_with_constraints(filtered_data, d, fc; init=p2v(d))
# end

# function fit_llh_with_constraints(data, f, fc; init = error("init!!"))
#     nll(p) = -sum((v>0) ? log(v) : -1e4 for v in f(data, p))
#     nll_constr(p) = nll(p) + fc(p)
#     return optimize_get_state_hessian(nll_constr, init, BFGS())
# end
