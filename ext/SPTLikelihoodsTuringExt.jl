"""
    SPTLikelihoodsTuringExt

Turing.jl support for the SPT-3G D1 T&E likelihood.

Loaded automatically when `SPTLikelihoods`, `Turing` and `Distributions` are all
present. It provides [`SPT3GD1Bandpowers`](@ref), a `Distribution` over the 1392
released bandpowers, so a model states the likelihood with `~`:

```julia
like.data.data_vector ~ SPT3GD1Bandpowers(like, model)
```

rather than injecting a number with `@addlogprob!`. `~` keeps the observation
visible to DynamicPPL, so conditioning, prior sampling and log-density
decompositions all behave normally.

## This one gives up nothing

Across these likelihood packages the `~` formulation has to trade something
away, and which thing depends on what the release stores. This package keeps
`cov_chol`, the Cholesky factor of the covariance itself, so:

  * `logpdf` is **normalized**. `logdet Σ` is a pass over the factor's diagonal,
    so the density is `loglikelihood(like, model) + gaussian_normalization(like)`
    at no cost. CamSpec has to take its normalization before `potri!` destroys
    the factor; JointCMB cannot afford one at all.
  * `rand` **works**, drawing `μ + L z` directly from the stored factor. The
    CamSpec and JointCMB distributions refuse, because recovering a factor from
    a stored precision matrix would mean an `n³` factorization hidden inside
    `rand`.

`logpdf` calls this package's own `chi2`, so reverse mode stays on the
registered `_fixed_covariance_solve` rule rather than re-deriving the solve.

## Priors

[`SPT3G_D1_PRIOR_DISTRIBUTIONS`](@ref) restates the Gaussian priors of the
default D1 candl YAML — the same 26 constraints `prior_penalty` applies, from
the same numbers. The test suite pins that the two agree up to a constant at
several parameter points, which is what keeps this from becoming a second
source of truth that can drift.

The remaining 17 parameters are unconstrained in that YAML. A Turing model
cannot leave them improper, so they carry deliberately wide, **unofficial**
defaults; `SPT3G_D1_UNCONSTRAINED_PARAMETERS` names them so a caller can see
exactly which ones are invented and replace them.
"""
module SPTLikelihoodsTuringExt

using SPTLikelihoods
using Turing
using Distributions
using LinearAlgebra
using Random

import SPTLikelihoods: chi2, predict, loglikelihood, gaussian_normalization

export SPT3GD1Bandpowers, spt3g_d1_model

"""
    SPT3GD1Bandpowers(like, model)

The SPT-3G D1 bandpower likelihood as a multivariate distribution: a normalized
Gaussian with the released fixed covariance and mean `model`.

The covariance is released data, never a parameter, so it is carried by
reference and its stored Cholesky factor is reused on every evaluation.
"""
struct SPT3GD1Bandpowers{L<:SPT3GD1Likelihood,M<:AbstractVector} <:
       ContinuousMultivariateDistribution
    like::L
    model::M

    function SPT3GD1Bandpowers(like::SPT3GD1Likelihood, model::AbstractVector)
        length(model) == length(like.data.data_vector) || throw(DimensionMismatch(
            "model must have $(length(like.data.data_vector)) entries, got $(length(model))",
        ))
        return new{typeof(like), typeof(model)}(like, model)
    end
end

Base.length(d::SPT3GD1Bandpowers) = length(d.like.data.data_vector)
Base.eltype(::Type{<:SPT3GD1Bandpowers{<:Any,M}}) where {M} = eltype(M)

# `chi2` measures from `like.data.data_vector`, so evaluating at an arbitrary
# `x` means shifting the mean by the same amount rather than calling it
# directly: chi2(like, model + (data - x)) is the quadratic form in (x - model).
function Distributions._logpdf(d::SPT3GD1Bandpowers, x::AbstractVector{<:Real})
    shifted = d.model .+ (d.like.data.data_vector .- x)
    return -chi2(d.like, shifted) / 2 + gaussian_normalization(d.like)
end

# The covariance factor is stored, so a draw is one triangular product. This is
# what makes prior predictive checks possible here.
function Distributions._rand!(rng::Random.AbstractRNG, d::SPT3GD1Bandpowers,
                              x::AbstractVector{<:Real})
    z = randn(rng, length(d))
    x .= d.model .+ d.like.data.cov_chol.L * z
    return x
end

"""
    SPT3G_D1_GAUSSIAN_PRIORS

The `(mean, sigma)` pairs of the default D1 candl YAML, for the 26 parameters it
constrains. Taken from the same numbers `prior_penalty` applies; the nine beam
modes are standard normal by construction.
"""
const SPT3G_D1_GAUSSIAN_PRIORS = (
    TT_CIBClustering_Alpha = (0.53, 0.1),
    EE_PolGalDust_Amp = (0.05, 0.022),
    EE_PolGalDust_Alpha = (-2.42, 0.04),
    EE_PolGalDust_Beta = (1.51, 0.04),
    TE_PolGalDust_Amp = (0.12, 0.051),
    TE_PolGalDust_Alpha = (-2.42, 0.04),
    TE_PolGalDust_Beta = (1.51, 0.04),
    TT_GalCirrus_Amp = (1.88, 0.96),
    TT_GalCirrus_Alpha = (-2.53, 0.05),
    TT_GalCirrus_Beta = (1.48, 0.02),
    TT_tSZ_Amp = (3.2279, 2.3764),
    TT_kSZ_Amp = (3.7287, 4.644),
    Kappa = (0.0, 0.00045),
    Tcal_ext150 = (1.0, 0.00360),
    beta_1 = (0.0, 1.0), beta_2 = (0.0, 1.0), beta_3 = (0.0, 1.0),
    beta_4 = (0.0, 1.0), beta_5 = (0.0, 1.0), beta_6 = (0.0, 1.0),
    beta_7 = (0.0, 1.0), beta_8 = (0.0, 1.0), beta_9 = (0.0, 1.0),
    T2P2_90 = (-0.006490, 0.001054),
    T2P2_150 = (-0.011941, 0.002103),
    T2P2_220 = (-0.022684, 0.006641),
)

"""
    SPT3G_D1_UNCONSTRAINED_PARAMETERS

The 17 parameters the default D1 candl YAML leaves unconstrained: the six
Poisson amplitudes, the three CIB amplitudes, the three polarized beam modes and
the five relative calibrations. Their entries in
[`SPT3G_D1_PRIOR_DISTRIBUTIONS`](@ref) are **not** official.
"""
const SPT3G_D1_UNCONSTRAINED_PARAMETERS = (
    :TT_Poisson_90x90, :TT_Poisson_90x150, :TT_Poisson_90x220,
    :TT_Poisson_150x150, :TT_Poisson_150x220, :TT_Poisson_220x220,
    :TT_CIB_150x150, :TT_CIB_150x220, :TT_CIB_220x220,
    :beta_pol_90, :beta_pol_150, :beta_pol_220,
    :Tcal_rel90, :Tcal_rel220, :Ecal_ext150, :Ecal_rel90, :Ecal_rel220,
)

# Wide, unofficial stand-ins for the parameters the YAML does not constrain.
# The amplitudes are positive quantities given ranges that comfortably contain
# the released posterior; the calibrations sit near unity.
const _UNCONSTRAINED_DEFAULTS = (
    TT_Poisson_90x90 = Uniform(0.0, 50.0),
    TT_Poisson_90x150 = Uniform(0.0, 50.0),
    TT_Poisson_90x220 = Uniform(0.0, 100.0),
    TT_Poisson_150x150 = Uniform(0.0, 50.0),
    TT_Poisson_150x220 = Uniform(0.0, 100.0),
    TT_Poisson_220x220 = Uniform(0.0, 200.0),
    TT_CIB_150x150 = Uniform(0.0, 50.0),
    TT_CIB_150x220 = Uniform(0.0, 100.0),
    TT_CIB_220x220 = Uniform(0.0, 200.0),
    beta_pol_90 = Normal(0.0, 1.0),
    beta_pol_150 = Normal(0.0, 1.0),
    beta_pol_220 = Normal(0.0, 1.0),
    Tcal_rel90 = Normal(1.0, 0.01),
    Tcal_rel220 = Normal(1.0, 0.01),
    Ecal_ext150 = Normal(1.0, 0.01),
    Ecal_rel90 = Normal(1.0, 0.01),
    Ecal_rel220 = Normal(1.0, 0.01),
)

const SPT3G_D1_PARAMETER_SYMBOLS = Symbol.(SPT3G_D1_PARAMETER_NAMES)

"""
    SPT3G_D1_PRIOR_DISTRIBUTIONS

One `Distribution` per entry of `SPT3G_D1_PARAMETER_NAMES`, in that order: the
official Gaussians of [`SPT3G_D1_GAUSSIAN_PRIORS`](@ref) where the YAML
constrains a parameter, and a wide unofficial default for each of
[`SPT3G_D1_UNCONSTRAINED_PARAMETERS`](@ref).
"""
const SPT3G_D1_PRIOR_DISTRIBUTIONS = NamedTuple{SPT3G_D1_PARAMETER_SYMBOLS}(
    map(SPT3G_D1_PARAMETER_SYMBOLS) do name
        gaussian = get(SPT3G_D1_GAUSSIAN_PRIORS, name, nothing)
        gaussian === nothing ? _UNCONSTRAINED_DEFAULTS[name] :
            Normal(gaussian[1], gaussian[2])
    end,
)

@eval @model function spt3g_d1_model(
    like::SPT3GD1Likelihood,
    foreground_model::SPT3GD1ForegroundModel,
    cmb::SPT3GD1CMBTheory;
    priors = SPT3G_D1_PRIOR_DISTRIBUTIONS,
)
    $(Expr(:block, [:($(name) ~ priors.$(name)) for name in SPT3G_D1_PARAMETER_SYMBOLS]...))
    parameters = SPT3GD1Parameters([$(SPT3G_D1_PARAMETER_SYMBOLS...)])

    model = predict(like, foreground_model, cmb, parameters)

    like.data.data_vector ~ SPT3GD1Bandpowers(like, model)

    return parameters
end

@doc """
    spt3g_d1_model(like, foreground_model, cmb; priors = SPT3G_D1_PRIOR_DISTRIBUTIONS)

Turing model for the SPT-3G D1 likelihood: the 43 nuisance parameters and then
the bandpower likelihood, all through `~`.

The parameter names and their order come from `SPT3G_D1_PARAMETER_NAMES`, so the
model cannot drift from the struct the package builds from that same vector.

The CMB spectra are held fixed; sample over cosmology by passing a different
`cmb` per evaluation or by wrapping this model in a larger one.
"""
spt3g_d1_model

end # module
