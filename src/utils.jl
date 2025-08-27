function compute_theory(DL_TT, DL_TE, DL_EE, κ,#1
    D_TT_90_90, D_TT_90_150, D_TT_90_220, D_TT_150_150, D_TT_150_220, D_TT_220_220,#6+1=7
    D_EE_90_90, D_EE_90_150, D_EE_90_220, D_EE_150_150, D_EE_150_220, D_EE_220_220,#6+7=13
    A_80_cirrus, α_cirrus, β_cirrus, A_80_cib, α_cib, β_cib, A_tSZ, ξ_tsz_CIB,#8+13=21-1=20 since α_cib_fixed
    A_kSZ, A_80_EE, α_EE, β_EE, A_80_TE, α_TE, β_TE,#7+20=27
    cal_T_90, cal_T_150, cal_T_220, cal_E_90, cal_E_150, cal_E_220, ℓs)#27+6=33

    ssl_TT = ssl_response(ℓs, κ, DL_TT)
    ssl_TE = ssl_response(ℓs, κ, DL_TE)
    ssl_EE = ssl_response(ℓs, κ, DL_EE)
    ab_coeff = -0.0004826
    ab_TT = aberration_response(ℓs, ab_coeff, DL_TT.+ssl_TT)
    ab_TE = aberration_response(ℓs, ab_coeff, DL_TE.+ssl_TE)
    ab_EE = aberration_response(ℓs, ab_coeff, DL_EE.+ssl_EE)
    ν_eff = effective_band_centres

    # Channel index order: 1→90 GHz, 2→150 GHz, 3→220 GHz
    pairs = ((1,1), (1,2), (1,3), (2,2), (2,3), (3,3))

    # Base spectra (identical across frequency pairs for a given type)
    base_TT = DL_TT .+ ssl_TT .+ ab_TT
    base_TE = DL_TE .+ ssl_TE .+ ab_TE
    base_EE = DL_EE .+ ssl_EE .+ ab_EE

    # Frequency lookups (non-mutating views or plain indexing are both fine)
    ν_gc_T = @view ν_eff[1, :]   # TT galactic-dust freqs
    ν_E    = @view ν_eff[2, :]   # E-channel freqs (EE/TE dust)
    ν_cib  = @view ν_eff[3, :]   # CIB freqs
    ν_tsz  = @view ν_eff[5, :]   # tSZ freqs

    # Shot-noise amplitudes (symmetric 3×3)
    D_TT = [D_TT_90_90   D_TT_90_150  D_TT_90_220;
            D_TT_90_150  D_TT_150_150 D_TT_150_220;
            D_TT_90_220  D_TT_150_220 D_TT_220_220]

    D_EE = [D_EE_90_90   D_EE_90_150  D_EE_90_220;
            D_EE_90_150  D_EE_150_150 D_EE_150_220;
            D_EE_90_220  D_EE_150_220 D_EE_220_220]

    # Calibration means (precompute with comprehensions; no mutation)
    cal_T = (cal_T_90, cal_T_150, cal_T_220)
    cal_E = (cal_E_90, cal_E_150, cal_E_220)

    cal_mean_TT = [cross_calibration_mean(cal_T[i], cal_T[j], cal_T[j], cal_T[i]) for i=1:3, j=1:3]
    cal_mean_TE = [cross_calibration_mean(cal_T[i], cal_E[j], cal_T[j], cal_E[i]) for i=1:3, j=1:3]
    cal_mean_EE = [cross_calibration_mean(cal_E[i], cal_E[j], cal_E[j], cal_E[i]) for i=1:3, j=1:3]

    # Foreground builders matching your existing signatures (pure, out-of-place)
    tt_fore(i,j) = TT_foregrounds(
        D_TT[i,j],
        A_80_cirrus, α_cirrus, β_cirrus, ν_gc_T[i], ν_gc_T[j],
        A_80_cib,    α_cib,    β_cib,    ν_cib[i],  ν_cib[j], 1.0, 1.0,
        A_tSZ,       ν_tsz[i], ν_tsz[j],
        ξ_tsz_CIB,   A_kSZ,    ℓs
    )

    ee_fore(i,j) = EE_foregrounds(
        D_EE[i,j],
        A_80_EE, α_EE, β_EE, ν_E[i], ν_E[j], ℓs
    )

    te_fore(i,j) = TE_foregrounds(
        A_80_TE, α_TE, β_TE, ν_E[i], ν_E[j], ℓs
    )

    # Final spectra per pair (pure, broadcasted arithmetic)
    tt_spec(i,j) = (base_TT .+ tt_fore(i,j)) ./ cal_mean_TT[i,j]
    te_spec(i,j) = (base_TE .+ te_fore(i,j)) ./ cal_mean_TE[i,j]
    ee_spec(i,j) = (base_EE .+ ee_fore(i,j)) ./ cal_mean_EE[i,j]

    # Assemble in the exact original ordering: TT_.., TE_.., EE_.. for each pair
    pure_theory = [f(i,j) for (i,j) in pairs for f in (tt_spec, te_spec, ee_spec)]

    model_matrix = stack(1:18) do i
       @view(window[:,:,i])' * pure_theory[i]
    end

    return model_matrix'
end

function slice_theory(model_matrix)
    residuals = bandpowers .- model_matrix
    vec_residuals = vcat([residuals[i, spec_bin_min[i]:spec_bin_max[i]] for i in 1:18]...)

    return vec_residuals
end

function compute_cov(model_matrix)
    dbs = vcat([model_matrix[i, spec_bin_min[i]:spec_bin_max[i]] for i in 1:18]...)
    Σ = cov .+ beam_cov .* kron(dbs, dbs')
    return Σ
end
