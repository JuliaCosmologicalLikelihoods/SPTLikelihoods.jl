module SPTLikelihoods

include("foregrounds.jl")

# Physical constants
T_CMB = 2.72548  # CMB temperature
h = 6.62606957e-34  # Planck's constant
kB = 1.3806488e-23  # Boltzmann constant
Ghz_Kelvin = h / kB * 1e9

end # module SPTLikelihoods
