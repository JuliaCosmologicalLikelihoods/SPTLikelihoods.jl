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

function GalacticDust(pow_at_80, α, β, ν1, ν2, SPT3G_windows_lmax, T_galdust, ν_0_galdust)

        # Grab ells helper (1-3200)
        ells = Array(1:SPT3G_windows_lmax)

        # Calculate and add galactic dust power
        Dl_galdust =  (ells ./ 80) .^ (α + 2.0) .* (pow_at_80 *
        DustFreqScaling(β, T_galdust, ν_0_galdust, ν1) *
        DustFreqScaling(β, T_galdust, ν_0_galdust, ν2))


        return Dl_galdust
end

function CIBClustering(pow_at_3000, α, β, ν1, ν2, z1, z2, SPT3G_windows_lmax, T_CIB,
    ν0_CIB)

        # Grab ells helper (1-3200)
        ells = Array(1:SPT3G_windows_lmax)

        # Calculate and add polarised galactic dust power
        Dl_cib_clustering =  (ells ./ 3000) .^ α .* (pow_at_3000 *
        DustFreqScaling(β, T_CIB, ν0_CIB, ν1) *
        DustFreqScaling(β, T_CIB, ν0_CIB, ν2) * sqrt(z1 * z2))

        return Dl_cib_clustering
end

function tSZCIBCorrelation(tSZ_template, ξ_tsz_CIB, tsz_pow_at_3000, CIB_pow_at_3000, α, β, z1, z2,
    CIB_ν1, CIB_ν2, tSZ_ν1, tSZ_ν2, SPT3G_windows_lmax, T_CIB, ν0_CIB, ν0_tSZ)

        # Calculate CIB components
        Dl_cib_clustering_11 = CIBClustering(
        CIB_pow_at_3000, α, β, CIB_ν1, CIB_ν1, z1, z1, SPT3G_windows_lmax, T_CIB, ν0_CIB)
        Dl_cib_clustering_22 = CIBClustering(
            CIB_pow_at_3000, α, β, CIB_ν2, CIB_ν2, z2, z2, SPT3G_windows_lmax, T_CIB, ν0_CIB)

        # Calculate the tSZ components
        Dl_tSZ_11 = tSZ(tSZ_template, tsz_pow_at_3000, tSZ_ν1, tSZ_ν1, ν0_tSZ)
        Dl_tSZ_22 = tSZ(tSZ_template, tsz_pow_at_3000, tSZ_ν2, tSZ_ν2, ν0_tSZ)

        # Calculate tSZ-CIB correlation
        # Sign defined such that a positive xi corresponds to a reduction at 150GHz
        #with np.errstate(invalid="ignore"):
        Dl_tSZ_CIB_corr = ( -1 * ξ_tsz_CIB
                .* (sqrt.(abs.(Dl_tSZ_11 .* Dl_cib_clustering_22)) .+
                   sqrt.(abs.(Dl_tSZ_22 .* Dl_cib_clustering_11))))

        return Dl_tSZ_CIB_corr
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
