struct OptimFitResults{T<:Optim.OptimizationResults, V<:Optim.AbstractOptimizerState} 
    optres::T
    state::V
    # forwarddiff_hessian_callback::Function
end

# the optimizer interface <: Optimizer should implement
#  - minimize(objective, init, optimizer)
#  - minimizer(minimizationResult)
#  - measurement(minimizationResult)
#  - minimum(minimizationResult)
struct BFGSApproxHesse end

function minimize(func, init, optimizer::BFGSApproxHesse;
        algorithm = BFGS(), autodiff = :forwarddiff, show_trace::Bool = true, kws...)
    obj = OnceDifferentiable(func, init; autodiff)
    # 
    options = Optim.Options(; show_trace, kws...)
    state = Optim.initial_state(algorithm, options, obj, init)
    #
    optres = optimize(obj, init, algorithm, options, state)
    #
    return OptimFitResults(optres, state)
end
#
# 
minimizer(fr::OptimFitResults) = minimizer(fr.optres)
minimum(fr::OptimFitResults) = minimum(fr.optres)
#
# 
invH(fr::OptimFitResults) = fr.state.invH
covmat(fr::OptimFitResults) = invH(fr)
errors(fr::OptimFitResults) = sqrt.(diag(fr.state.invH))
#
hessian_callback(p) = ForwardDiff.hessian(func, p)
invexacthessian(fr::OptimFitResults, p = minimizer(fr)) = inv(fr.forwarddiff_hessian_callback(p))
#
# 
function measurements(fr, exacthessian::Bool=false)
    mv = minimizer(fr)
    δv = !(exacthessian) ? errors(fr) : sqrt.(diag(invexacthessian(fr)))
    return [±(m,δ) for (m,δ) in zip(mv, δv)]
end
