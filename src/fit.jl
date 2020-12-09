#
llh(data, f; weights = fill(1.0, length(data))) = -sum((v>0) ? w*log(v) : -1e4 for (w,v) in zip(weights, f(data)))

function fit_llh(data, f; init_pars = error("init_pars!!"), weights = fill(1.0, length(data)))
    llh(p) = -sum((v>0) ? w*log(v) : -1e4 for (w,v) in zip(weights, f(data, p)))
    llh(init_pars) # test
    pfr = Optim.optimize(llh, init_pars, BFGS(),
                   Optim.Options(show_trace = true); autodiff = :forwarddiff)
    return pfr
end

function fit_llh(data, d::T where T<:AdvancedFunction)
    filtered_data = filter(x->inrange(x, lims(d)), data)
    return fit_llh(filtered_data, d; init_pars=p2v(d))
end
