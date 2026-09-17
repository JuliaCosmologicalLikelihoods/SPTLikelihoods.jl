function _gaussian_penalty(value, mean, sigma)
    return ((value - mean) / sigma)^2 / 2
end

"""
    prior_penalty(parameters; tau=nothing)

Return the unnormalized positive Gaussian-prior penalty in the default D1
candl YAML. The data likelihood remains prior-free. `tau` is optional because
standard joint analyses clear candl's internal tau prior and impose it globally.
"""
function prior_penalty(parameters::SPT3GD1Parameters; tau=nothing)
    penalty = zero(parameters.kappa)
    penalty += _gaussian_penalty(parameters.cib_alpha, 0.53, 0.1)
    penalty += _gaussian_penalty(parameters.dust_ee.amplitude, 0.05, 0.022)
    penalty += _gaussian_penalty(parameters.dust_ee.alpha, -2.42, 0.04)
    penalty += _gaussian_penalty(parameters.dust_ee.beta, 1.51, 0.04)
    penalty += _gaussian_penalty(parameters.dust_te.amplitude, 0.12, 0.051)
    penalty += _gaussian_penalty(parameters.dust_te.alpha, -2.42, 0.04)
    penalty += _gaussian_penalty(parameters.dust_te.beta, 1.51, 0.04)
    penalty += _gaussian_penalty(parameters.dust_tt.amplitude, 1.88, 0.96)
    penalty += _gaussian_penalty(parameters.dust_tt.alpha, -2.53, 0.05)
    penalty += _gaussian_penalty(parameters.dust_tt.beta, 1.48, 0.02)
    penalty += _gaussian_penalty(parameters.tsz, 3.2279, 2.3764)
    penalty += _gaussian_penalty(parameters.ksz, 3.7287, 4.644)
    penalty += _gaussian_penalty(parameters.kappa, 0.0, 0.00045)
    penalty += _gaussian_penalty(parameters.tcal_ext150, 1.0, 0.00360)
    for beta in parameters.beam_modes
        penalty += beta^2 / 2
    end
    penalty += _gaussian_penalty(parameters.t2p[1], -0.006490, 0.001054)
    penalty += _gaussian_penalty(parameters.t2p[2], -0.011941, 0.002103)
    penalty += _gaussian_penalty(parameters.t2p[3], -0.022684, 0.006641)
    if !isnothing(tau)
        penalty += _gaussian_penalty(tau, 0.051, 0.006)
    end
    return penalty
end

prior_chi2(parameters::SPT3GD1Parameters; tau=nothing) = 2 * prior_penalty(parameters; tau=tau)
logprior(parameters::SPT3GD1Parameters; tau=nothing) = -prior_penalty(parameters; tau=tau)

function logposterior(
    like::SPT3GD1Likelihood,
    model_bandpowers::AbstractVector,
    parameters::SPT3GD1Parameters;
    tau=nothing,
)
    return loglikelihood(like, model_bandpowers) + logprior(parameters; tau=tau)
end
