Base.@kwdef struct SPT3G_2018_TTTEEE_Foregrounds
    ℓ_min::Integer=1
    ℓ_max::Integer=3200
    galdust_ν0=150
    galdust_T=19.6
    CIB_ν0 = 150.
    CIB_T = 25.
    tSZ_ν0 = 143
end

#Planck function normalised to 1 at ν0
function _Bnu(ν, ν0, T)
    return (ν / ν0)^3 * (exp(Ghz_Kelvin * ν0 / T) - 1.0) / (exp(Ghz_Kelvin * ν / T) - 1.0)
end

# Derivative of Planck function normalised to 1 at ν0
function _dBdT(ν, ν0, T)
    x0 = Ghz_Kelvin * ν0 / T
    x = Ghz_Kelvin * ν / T

    dBdT0 = x0^4 * exp(x0) / (exp(x0) - 1)^2
    dBdT = x^4 * exp(x) / (exp(x) - 1)^2

    return dBdT / dBdT0
end

function _dust_f_scaling(β, Tdust, ν0, ν_eff)
    return (ν_eff / ν0) ^ β * _Bnu(ν_eff, ν0, Tdust) / _dBdT(ν_eff, ν0, T_CMB)
end

function _galactic_dust(A_80, α, β, ν1, ν2, SPT3G_windows_lmax)

    ells = Array{Float64}(1:SPT3G_windows_lmax)
    #TODO optimization opportunity:  the α exponent has a small dynamical range in the
    # chains. Maybe an interpolator? Also, probably not doing the best thing even without it

    # Calculate and add galactic dust power
    Dl_galdust =  (ells ./ 80) .^ (α + 2.0) .* (A_80 *
    _dust_f_scaling(β, galdust_T, galdust_ν0, ν1) *
    _dust_f_scaling(β, galdust_T, galdust_ν0, ν2))

    return Dl_galdust
end

function _CIB_clustering(pow_at_3000, α, β, ν1, ν2, z1, z2, SPT3G_windows_lmax)

        ells = Array{Float64}(1:SPT3G_windows_lmax)

        # Calculate and add polarised galactic dust power
        Dl_cib_clustering =  (ells ./ 3000) .^ α .* (pow_at_3000 *
        _dust_f_scaling(β, CIB_T, CIB_ν0, ν1) *
        _dust_f_scaling(β, CIB_T, CIB_ν0, ν2) * sqrt(z1 * z2))

        return Dl_cib_clustering
end

function _tSZ_CIB_correlation(ξ_tsz_CIB, tsz_pow_at_3000, CIB_pow_at_3000, α, β,
    z1, z2, CIB_ν1, CIB_ν2, tSZ_ν1, tSZ_ν2, SPT3G_windows_lmax)

        # Calculate CIB components
        Dl_cib_clustering_11 = _CIB_clustering(
        CIB_pow_at_3000, α, β, CIB_ν1, CIB_ν1, z1, z1, SPT3G_windows_lmax)
        Dl_cib_clustering_22 = _CIB_clustering(
            CIB_pow_at_3000, α, β, CIB_ν2, CIB_ν2, z2, z2, SPT3G_windows_lmax)

        # Calculate the tSZ components
        Dl_tSZ_11 = _tSZ(tsz_pow_at_3000, tSZ_ν1, tSZ_ν1)
        Dl_tSZ_22 = _tSZ(tsz_pow_at_3000, tSZ_ν2, tSZ_ν2)

        # Calculate tSZ-CIB correlation
        # Sign defined such that a positive xi corresponds to a reduction at 150GHz
        #TODO: maybe use NaNMath.jl to deal with possible NaNs?
        """Dl_tSZ_CIB_corr = ( -1 * ξ_tsz_CIB
                .* (sqrt.(abs.(Dl_tSZ_11 .* Dl_cib_clustering_22)) .+
                   sqrt.(abs.(Dl_tSZ_22 .* Dl_cib_clustering_11))))"""

        return ( -1 * ξ_tsz_CIB .* (sqrt.(Dl_tSZ_11 .* Dl_cib_clustering_22) .+
                                    sqrt.(Dl_tSZ_22 .* Dl_cib_clustering_11)))
end

function _tSZ_f_scaling(ν, ν0, T)
    x0 = Ghz_Kelvin * ν0 / T
    x = Ghz_Kelvin * ν / T

    tSZfac0 = x0 * (exp(x0) + 1) / (exp(x0) - 1) - 4
    tSZfac = x * (exp(x) + 1) / (exp(x) - 1) - 4
    #negative for ν = 219

    return tSZfac / tSZfac0
end


function _tSZ(A_tSZ, ν1, ν2)

    # Calculate tSZ power
    Dl_tSZ = tSZ_template .* (A_tSZ *
    _tSZ_f_scaling(ν1, tSZ_ν0, T_CMB) *
    _tSZ_f_scaling(ν2, tSZ_ν0, T_CMB))
    # Frequency scaling

    return Dl_tSZ
end

function _kSZ(pow_at_3000)
    return pow_at_3000 .* kSZ_template
end

function _getCℓ_derivative(SPT3G_windows_lmax, Dℓ_theory)
    ells = Array(1:SPT3G_windows_lmax)

    Cℓ_derivative = Dℓ_theory * 2 * π ./ (ells .* (ells .+ 1))
    Cℓ_derivative[2:end-1] .= 0.5 .* (Cℓ_derivative[3:end]-Cℓ_derivative[1:end-2])
    Cℓ_derivative[1] = Cℓ_derivative[2]
    Cℓ_derivative[end] = Cℓ_derivative[end-1]

    return Cℓ_derivative
end

function _supersamplelensing(SPT3G_windows_lmax, κ, Dℓ_theory)
    ells = Array(1:SPT3G_windows_lmax)

    Cℓ_derivative = _getCℓ_derivative(SPT3G_windows_lmax, Dℓ_theory)
    ssl_correction = ells .* Cℓ_derivative  .* ells .* (ells .+ 1) ./ (2π)
    #maybe better to have a get_Dl_derivative?
    ssl_correction .+= 2 .* Dℓ_theory
    ssl_correction .*= (-κ)

    #TODO: pay attention: this is named apply, but you are not applying the correction
    # you are just evaluating it!

    return ssl_correction
end

function _abberation_correction(SPT3G_windows_lmax, ab_coeff, Dℓ_theory)
    ells = Array(1:SPT3G_windows_lmax)

    Cℓ_derivative = _getCℓ_derivative(SPT3G_windows_lmax, Dℓ_theory)
    aberration_correction = -ab_coeff .* Cℓ_derivative .* ells
    aberration_correction .*= ells .* (ells .+ 1) ./ (2π)

    #TODO: pay attention: this is named apply, but you are not applying the correction
    # you are just evaluating it!

    return aberration_correction
end

function _calibration(cal1, cal2, cal3, cal4)
    return 0.5 * (cal1 * cal2 + cal3 * cal4)
end

function _poisson_power(SPT3G_windows_lmax, pow_at_3000)
    ells = Array(1:SPT3G_windows_lmax)

    return ells .* ells .* (pow_at_3000 / 3000 ^2)
end

### implementation proposal
### write a function for each of the spectra type
### in your likelihood, you are gonna call them several times. Beautiful? Not. But it should work

function TT_foregrounds(D_TT_ν1_ν2, A_80_cirrus, α_cirrus, β_cirrus,
    ν1_gc, ν2_gc, A_80_cib, α_cib, β_cib, ν1_cib, ν2_cib, z1, z2, A_tSZ, ν1_tSZ, ν2_tSZ,
    ξ_tsz_CIB, A_kSZ, SPT3G_windows_lmax)###TT_params
    #α_cib fixed to 0.8, we will do it at the likelihood level
    #z1, z2 fixed to 1.
    pp = _poisson_power(SPT3G_windows_lmax, D_TT_ν1_ν2)
    gd = _galactic_dust(A_80_cirrus, α_cirrus, β_cirrus, ν1_gc, ν2_gc, SPT3G_windows_lmax)
    cc = _CIB_clustering(A_80_cib, α_cib, β_cib, ν1_cib, ν2_cib, z1, z2, SPT3G_windows_lmax)
    tsz = _tSZ(A_tSZ, ν1_tSZ, ν2_tSZ)
    tszcib = _tSZ_CIB_correlation(ξ_tsz_CIB, A_tSZ, A_80_cib, α_cirrus, β_cirrus,
    z1, z2, ν1_cib, ν2_cib, ν1_tSZ, ν2_tSZ, SPT3G_windows_lmax)
    ksz = _kSZ(A_kSZ)
    return pp.+gd.+cc.+tsz.+tszcib.+ksz
end

function TE_foregrounds(A_80_TE, α_TE, β_TE, ν1_te, ν2_te, SPT3G_windows_lmax)###TE_params
    return _galactic_dust(A_80_TE, α_TE, β_TE, ν1_te, ν2_te, SPT3G_windows_lmax)
end

function EE_foregrounds(D_EE_ν1_ν2, A_80_EE, α_EE, β_EE, ν1_ee, ν2_ee, SPT3G_windows_lmax)###EE_params
    pp = _poisson_power(SPT3G_windows_lmax, D_EE_ν1_ν2)
    gd = _galactic_dust(A_80_EE, α_EE, β_EE, ν1_ee, ν2_ee, SPT3G_windows_lmax)
    return pp+gd
end
