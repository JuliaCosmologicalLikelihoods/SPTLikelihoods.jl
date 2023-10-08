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

end


include("foregrounds.jl")

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

end # module SPTLikelihoods
