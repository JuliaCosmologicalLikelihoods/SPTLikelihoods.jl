module SPTLikelihoods
import Base.@kwdef
using Artifacts
using LinearAlgebra
using CMBForegrounds: CMBForegrounds, window_convolution, correlation_power, GeometricMeanCorrelation,
    PowerLawShape, angular_power
using NPZ
using DelimitedFiles

export SPT3GD1Data, load_spt3g_d1_data, SPT3GD1Likelihood, bin_theory, chi2, loglikelihood
export SPT3GD1CMBTheory, SPT3GD1ForegroundModel, SPT3GD1Parameters, SPT3G_D1_PARAMETER_NAMES, SPT3G_D1_FIDUCIAL_PARAMETERS, foreground_stages
export instrument_stages, predict
export prior_penalty, prior_chi2, logprior, logposterior
export compute_theory, slice_theory, compute_cov

const tSZ_template = Vector{Float64}(undef, 3200)
const kSZ_template = Vector{Float64}(undef, 3200)
const window = Array{Float64,3}(undef, 3200, 44, 18)
const bandpowers = Matrix{Float64}(undef, 18, 44)
const beam_cov = Matrix{Float64}(undef, 728, 728)
const cal_cov = Matrix{Float64}(undef, 6, 6)
const effective_band_centres = Matrix{Float64}(undef, 5, 3)
const fid_cov = Matrix{Float64}(undef, 768, 768)
const cov = Matrix{Float64}(undef, 728, 728)

const _legacy_2018_loaded = Ref(false)

function _ensure_legacy_2018_loaded!()
    _legacy_2018_loaded[] && return
    tSZ_template .= npzread(joinpath(artifact"SPT3G_data",
        "tSZ_Dl_shaw10_153ghz_norm1.npy"))[1:3200, 2]
    kSZ_template .= npzread(joinpath(artifact"SPT3G_data",
        "kSZ_Dl_CSF_incl_patchy_norm1.npy"))[1:3200, 2]

    window .= permutedims(npzread(joinpath(artifact"SPT3G_data", "windows.npy")), (3, 1, 2))

    bandpowers .= npzread(joinpath(artifact"SPT3G_data",
        "SPT3G_2018_TTTEEE_bandpowers.npy"))
    beam_cov .= npzread(joinpath(artifact"SPT3G_data",
        "SPT3G_2018_TTTEEE_beam_covariance.npy"))
    cal_cov .= npzread(joinpath(artifact"SPT3G_data",
        "SPT3G_2018_TTTEEE_cal_covariance.npy"))
    effective_band_centres .= npzread(joinpath(artifact"SPT3G_data",
        "SPT3G_2018_TTTEEE_effective_band_centres.npy"))
    fid_cov .= npzread(joinpath(artifact"SPT3G_data",
        "SPT3G_2018_TTTEEE_fiducial_covariance.npy"))
    cov .= npzread(joinpath(artifact"SPT3G_data",
        "bp_cov_posdef.npy"))
    _legacy_2018_loaded[] = true
    return
end


include("foregrounds.jl")
include("utils.jl")
include("spt3g_d1_data.jl")
include("spt3g_d1_likelihood.jl")
include("spt3g_d1_foregrounds.jl")
include("spt3g_d1_instrument.jl")
include("spt3g_d1_priors.jl")

# Physical constants
const T_CMB = 2.72548  # CMB temperature
const h = 6.62606957e-34  # Planck's constant
const kB = 1.3806488e-23  # Boltzmann constant
const Ghz_Kelvin = h / kB * 1e9
const galdust_T = 19.6
const galdust_ν0 = 150.0
const CIB_ν0 = 150.0
const CIB_T = 25.0
const tSZ_ν0 = 143.0
const spec_bin_min = [10, 1, 1, 10, 1, 1, 10, 1, 1, 10, 1, 1, 15, 1, 1, 15, 1,
    1]
const spec_bin_max = [44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44,
    44]
const bin_min = 1
const bin_max = 44
const windows_lmin = 1
const windows_lmax = 3200


end # module SPTLikelihoods
