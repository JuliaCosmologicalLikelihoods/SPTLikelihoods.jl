#Planck distribution
function Bnu(ν, ν0, T)
    return (ν / ν0)^3 * (exp(Ghz_Kelvin * ν0 / T) - 1.0) / (exp(Ghz_Kelvin * ν / T) - 1.0)
end


# Derivative of Planck function normalised to 1 at nu0
function dBdT(ν, ν0, T)
    x0 = Ghz_Kelvin * ν0 / T
    x = Ghz_Kelvin * ν / T

    dBdT0 = x0^4 * exp(x0) / (exp(x0) - 1)^2
    dBdT = x^4 * exp(x) / (exp(x) - 1)^2

    return dBdT / dBdT0
end

function DustFreqScaling(β, Tdust, ν0, ν_eff)
    fdust = (ν_eff / ν0) ^ β * Bnu(ν_eff, ν0, Tdust) / dBdT(ν_eff, ν0, T_CMB)
    return fdust
end

function tSZFrequencyScaling(ν, ν0, T)
    x0 = Ghz_Kelvin * ν0 / T
    x = Ghz_Kelvin * ν / T

    tSZfac0 = x0 * (exp(x0) + 1) / (exp(x0) - 1) - 4
    tSZfac = x * (exp(x) + 1) / (exp(x) - 1) - 4

    return tSZfac / tSZfac0
end

function tSZ(tSZ_template, pow_at_3000, ν1, ν2, ν0_tSZ)

    # Calculate tSZ power
    Dl_tSZ = tSZ_template .* (pow_at_3000 *
    tSZFrequencyScaling(ν1, ν0_tSZ, T_CMB) *
    tSZFrequencyScaling(ν2, ν0_tSZ, T_CMB))
    # Frequency scaling

    return Dl_tSZ
end
