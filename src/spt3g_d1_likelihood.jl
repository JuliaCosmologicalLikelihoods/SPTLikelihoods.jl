"""
    SPT3GD1Likelihood(data)

Gaussian fixed-covariance likelihood for the full SPT-3G D1 T&E release.
Foreground and instrument transformations are added separately; this initial
boundary accepts the final release-ordered unbinned Dℓ vector.
"""
struct SPT3GD1Likelihood
    data::SPT3GD1Data
end

SPT3GD1Likelihood(data_dir::AbstractString) = SPT3GD1Likelihood(load_spt3g_d1_data(data_dir))

function SPT3GD1Likelihood()
    artifact_dir = _get_d1_artifact_dir()
    if isnothing(artifact_dir) || !isdir(artifact_dir)
        throw(ArgumentError("SPT-3G D1 artifact \"SPT3G_D1_TnE_v0\" is not bound or installed in Artifacts.toml. Please provide the artifact path via SPT3GD1Likelihood(data_dir)."))
    end
    return SPT3GD1Likelihood(load_spt3g_d1_data(artifact_dir; require_metadata=true))
end

@inline _fixed_covariance_solve(cov_chol::Cholesky, residual::AbstractVector) = cov_chol \ residual

"""
    bin_theory(like, unbinned_Dls)

Apply the 21 released window matrices to a spectrum-major unbinned Dℓ vector.
Each block must have the D1 theory-grid length and follow
`like.data.spectrum_order` exactly.
"""
function bin_theory(like::SPT3GD1Likelihood, unbinned_Dls::AbstractVector)
    data = like.data
    n_ell = length(data.ells)
    expected_length = length(data.windows) * n_ell
    length(unbinned_Dls) == expected_length ||
        throw(DimensionMismatch("unbinned Dℓ vector must have length $expected_length"))

    bandpowers = map(enumerate(data.windows)) do (index, window)
        block = (index - 1) * n_ell + 1:index * n_ell
        window_convolution(window, @view unbinned_Dls[block])
    end
    return reduce(vcat, bandpowers)
end

"""
    chi2(like, model_bandpowers)

Return `(data - model)' * covariance^-1 * (data - model)` using the fixed
lower Cholesky factor stored in `like`.
"""
function chi2(like::SPT3GD1Likelihood, model_bandpowers::AbstractVector)
    data = like.data
    length(model_bandpowers) == length(data.data_vector) ||
        throw(DimensionMismatch("model bandpowers must have length $(length(data.data_vector))"))
    residual = data.data_vector .- model_bandpowers
    solved = _fixed_covariance_solve(data.cov_chol, residual)
    return dot(residual, solved)
end

chi2(like::SPT3GD1Likelihood, unbinned_Dls::AbstractVector, ::Val{:unbinned}) =
    chi2(like, bin_theory(like, unbinned_Dls))

"""
    loglikelihood(like, model_bandpowers)

Return the unnormalized data-only log likelihood `-chi2/2`. The parameter
independent Gaussian normalization and external cosmological priors are omitted.
"""
loglikelihood(like::SPT3GD1Likelihood, model_bandpowers::AbstractVector) =
    -chi2(like, model_bandpowers) / 2

loglikelihood(like::SPT3GD1Likelihood, unbinned_Dls::AbstractVector, ::Val{:unbinned}) =
    -chi2(like, unbinned_Dls, Val(:unbinned)) / 2
