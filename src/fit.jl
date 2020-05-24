# 
function fit_llh(data, f; init_pars = error("init_pars!!"), weights = fill(1.0, length(data)))
    llh(p) = -sum((v>0) ? w*log(v) : -1e4 for (w,v) in zip(weights, f(data, p)))
    llh(init_pars) # test
    pfr = Optim.minimizer(Optim.optimize(llh, init_pars, BFGS(),
                   Optim.Options(show_trace = true); autodiff = :forwarddiff))
    return pfr
end
