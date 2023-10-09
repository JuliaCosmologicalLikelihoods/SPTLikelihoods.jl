module SPTLikelihoods
import Base.@kwdef
using Artifacts
using NPZ

function __init__()
    #check : lmax=3200?

    global tSZ_template = npzread(joinpath(artifact"SPT3G_data",
                                           "tSZ_Dl_shaw10_153ghz_norm1.npy"))[1:3200, 2]
    global kSZ_template = npzread(joinpath(artifact"SPT3G_data",
                                           "kSZ_Dl_CSF_incl_patchy_norm1.npy"))[1:3200, 2]
    global window = zeros(44,18,3200)
    for i in 1:44
        window[i,:,:] = npzread(joinpath(artifact"SPT3G_data", "window_"*string(i)*".npy"))
    end

    global new_window = npzread(joinpath(artifact"SPT3G_data", "windows.npy"))

    global bandpowers = npzread(joinpath(artifact"SPT3G_data",
                                           "SPT3G_2018_TTTEEE_bandpowers.npy"))
    global beam_cov = npzread(joinpath(artifact"SPT3G_data",
                                           "SPT3G_2018_TTTEEE_beam_covariance.npy"))
    global cal_cov = npzread(joinpath(artifact"SPT3G_data",
                                           "SPT3G_2018_TTTEEE_cal_covariance.npy"))
    global effective_band_centres = npzread(joinpath(artifact"SPT3G_data",
                                           "SPT3G_2018_TTTEEE_effective_band_centres.npy"))
    global fid_cov = npzread(joinpath(artifact"SPT3G_data",
                                           "SPT3G_2018_TTTEEE_fiducial_covariance.npy"))
    global cov = npzread(joinpath(artifact"SPT3G_data",
                                           "bp_cov_posdef.npy"))

end


include("foregrounds.jl")
include("utils.jl")

# Physical constants
const T_CMB = 2.72548  # CMB temperature
const h = 6.62606957e-34  # Planck's constant
const kB = 1.3806488e-23  # Boltzmann constant
const Ghz_Kelvin = h / kB * 1e9
const galdust_T = 19.6
const galdust_ν0 = 150.
const CIB_ν0 = 150.
const CIB_T = 25.0
const tSZ_ν0 = 143.
const spec_bin_min = [10,  1,  1, 10,  1,  1, 10,  1,  1, 10,  1,  1, 15,  1,  1, 15,  1,
                      1]
const spec_bin_max = [44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44, 44,
                      44]
const bin_min = 1
const bin_max = 44
const windows_lmin = 1
const windows_lmax = 3200


end # module SPTLikelihoods
