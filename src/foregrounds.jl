dust_tt_power_law(ℓs, A_80, α, β, ν1, ν2) = CMBForegrounds.dust_tt_power_law(ℓs, A_80, α, β, ν1, ν2, galdust_T, galdust_ν0; ℓ_pivot=80, T_CMB=T_CMB)

cib_clustered_power(ℓs, pow_at_3000, α, β, ν1, ν2, z1, z2) = CMBForegrounds.cib_clustered_power(ℓs, pow_at_3000, α, β, ν1, ν2, z1, z2, CIB_T, CIB_ν0)

tsz_cib_cross_power(ℓs, ξ_tsz_CIB, tsz_pow_at_3000, CIB_pow_at_3000, α, β,
    z1, z2, CIB_ν1, CIB_ν2, tSZ_ν1, tSZ_ν2) =
    CMBForegrounds.tsz_cib_cross_power(ℓs, ξ_tsz_CIB, tsz_pow_at_3000, CIB_pow_at_3000, α, β,
        z1, z2, CIB_ν1, CIB_ν2, tSZ_ν1, tSZ_ν2, tSZ_template, tSZ_ν0, CIB_T, CIB_ν0; ℓ_pivot=3000, T_CMB=T_CMB
    )

tsz_cross_power(A_tSZ, ν1, ν2) = CMBForegrounds.tsz_cross_power(tSZ_template, A_tSZ, ν1, ν2, tSZ_ν0)

ksz_template_scaled(pow_at_3000) = CMBForegrounds.ksz_template_scaled(kSZ_template, pow_at_3000)

ssl_response(ls, κ, Dl) = CMBForegrounds.ssl_response(ls, κ, Dl)

aberration_response(ℓs, ab_coeff, Dℓ_theory) = CMBForegrounds.aberration_response(ℓs, ab_coeff, Dℓ_theory)

cross_calibration_mean(cal1, cal2, cal3, cal4) = CMBForegrounds.cross_calibration_mean(cal1, cal2, cal3, cal4)

shot_noise_power(ℓs, pow_at_3000) = CMBForegrounds.shot_noise_power(ℓs, pow_at_3000; ℓ0=3000)

function TT_foregrounds(D_TT_ν1_ν2, A_80_cirrus, α_cirrus, β_cirrus,
    ν1_gc, ν2_gc, A_80_cib, α_cib, β_cib, ν1_cib, ν2_cib, z1, z2, A_tSZ, ν1_tSZ, ν2_tSZ,
    ξ_tsz_CIB, A_kSZ, ℓs)###TT_params
    #α_cib fixed to 0.8, we will do it at the likelihood level
    #z1, z2 fixed to 1.
    pp = shot_noise_power(ℓs, D_TT_ν1_ν2)
    gd = dust_tt_power_law(ℓs, A_80_cirrus, α_cirrus, β_cirrus, ν1_gc, ν2_gc)
    cc = cib_clustered_power(ℓs, A_80_cib, α_cib, β_cib, ν1_cib, ν2_cib, z1, z2)
    tsz = tsz_cross_power(A_tSZ, ν1_tSZ, ν2_tSZ)
    tszcib = tsz_cib_cross_power(ℓs, ξ_tsz_CIB, A_tSZ, A_80_cib, α_cib, β_cib,
        z1, z2, ν1_cib, ν2_cib, ν1_tSZ, ν2_tSZ)
    ksz = ksz_template_scaled(A_kSZ)
    return pp .+ gd .+ cc .+ tsz .+ tszcib .+ ksz
end

function TE_foregrounds(A_80_TE, α_TE, β_TE, ν1_te, ν2_te, ℓs)###TE_params
    return dust_tt_power_law(ℓs, A_80_TE, α_TE, β_TE, ν1_te, ν2_te)
end

function EE_foregrounds(D_EE_ν1_ν2, A_80_EE, α_EE, β_EE, ν1_ee, ν2_ee, ℓs)###EE_params
    pp = shot_noise_power(ℓs, D_EE_ν1_ν2)
    gd = dust_tt_power_law(ℓs, A_80_EE, α_EE, β_EE, ν1_ee, ν2_ee)
    return pp + gd
end
