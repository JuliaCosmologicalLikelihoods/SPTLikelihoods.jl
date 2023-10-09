function compute_theory(DL_TT, DL_TE, DL_EE, κ,#1
    D_TT_90_90, D_TT_90_150, D_TT_90_220, D_TT_150_150, D_TT_150_220, D_TT_220_220,#6+1=7
    D_EE_90_90, D_EE_90_150, D_EE_90_220, D_EE_150_150, D_EE_150_220, D_EE_220_220,#6+7=13
    A_80_cirrus, α_cirrus, β_cirrus, A_80_cib, α_cib, β_cib, A_tSZ, ξ_tsz_CIB,#8+13=21-1=20 since α_cib_fixed
    A_kSZ, A_80_EE, α_EE, β_EE, A_80_TE, α_TE, β_TE,#7+20=27
    cal_T_90, cal_T_150, cal_T_220, cal_E_90, cal_E_150, cal_E_220, SPT3G_windows_lmax)#27+6=33

    ssl_TT = _supersamplelensing(SPT3G_windows_lmax, κ, DL_TT)
    ssl_TE = _supersamplelensing(SPT3G_windows_lmax, κ, DL_TE)
    ssl_EE = _supersamplelensing(SPT3G_windows_lmax, κ, DL_EE)
    ab_coeff = -0.0004826
    ab_TT = _abberation_correction(SPT3G_windows_lmax, ab_coeff, DL_TT)
    ab_TE = _abberation_correction(SPT3G_windows_lmax, ab_coeff, DL_TE)
    ab_EE = _abberation_correction(SPT3G_windows_lmax, ab_coeff, DL_EE)
    ν_eff = effective_band_centres

    TT_90_90   = (DL_TT.+ssl_TT.+ab_TT.+TT_foregrounds(D_TT_90_90, A_80_cirrus, α_cirrus, β_cirrus,
    ν_eff[1,1], ν_eff[1,1], A_80_cib, α_cib, β_cib, ν_eff[3,1], ν_eff[3,1], 1., 1., A_tSZ,
    ν_eff[5,1], ν_eff[5,1], ξ_tsz_CIB, A_kSZ, SPT3G_windows_lmax)) ./
    _calibration(cal_T_90, cal_T_90, cal_T_90, cal_T_90)

    TT_90_150  = (DL_TT.+ssl_TT.+ab_TT.+TT_foregrounds(D_TT_90_150, A_80_cirrus, α_cirrus, β_cirrus,
    ν_eff[1,1], ν_eff[1,2], A_80_cib, α_cib, β_cib, ν_eff[3,1], ν_eff[3,2], 1., 1., A_tSZ,
    ν_eff[5,1], ν_eff[5,2], ξ_tsz_CIB, A_kSZ, SPT3G_windows_lmax)) ./
    _calibration(cal_T_90, cal_T_150, cal_T_90, cal_T_150)

    TT_90_220  = (DL_TT.+ssl_TT.+ab_TT.+TT_foregrounds(D_TT_90_220, A_80_cirrus, α_cirrus, β_cirrus,
    ν_eff[1,1], ν_eff[1,3], A_80_cib, α_cib, β_cib, ν_eff[3,1], ν_eff[3,3], 1., 1., A_tSZ,
    ν_eff[5,1], ν_eff[5,3], ξ_tsz_CIB, A_kSZ, SPT3G_windows_lmax)) ./
    _calibration(cal_T_90, cal_T_220, cal_T_90, cal_T_220)

    TT_150_150 = (DL_TT.+ssl_TT.+ab_TT.+TT_foregrounds(D_TT_150_150, A_80_cirrus, α_cirrus, β_cirrus,
    ν_eff[1,2], ν_eff[1,2], A_80_cib, α_cib, β_cib, ν_eff[3,2], ν_eff[3,2], 1., 1., A_tSZ,
    ν_eff[5,2], ν_eff[5,2], ξ_tsz_CIB, A_kSZ, SPT3G_windows_lmax)) ./
    _calibration(cal_T_150, cal_T_150, cal_T_150, cal_T_150)

    TT_150_220 = (DL_TT.+ssl_TT.+ab_TT.+TT_foregrounds(D_TT_150_220, A_80_cirrus, α_cirrus, β_cirrus,
    ν_eff[1,2], ν_eff[1,3], A_80_cib, α_cib, β_cib, ν_eff[3,2], ν_eff[3,3], 1., 1., A_tSZ,
    ν_eff[5,2], ν_eff[5,3], ξ_tsz_CIB, A_kSZ, SPT3G_windows_lmax)) ./
    _calibration(cal_T_150, cal_T_150, cal_T_220, cal_T_220)

    TT_220_220 = (DL_TT.+ssl_TT.+ab_TT.+TT_foregrounds(D_TT_220_220, A_80_cirrus, α_cirrus, β_cirrus,
    ν_eff[1,3], ν_eff[1,3], A_80_cib, α_cib, β_cib, ν_eff[3,3], ν_eff[3,3], 1., 1., A_tSZ,
    ν_eff[5,3], ν_eff[5,3], ξ_tsz_CIB, A_kSZ, SPT3G_windows_lmax)) ./
    _calibration(cal_T_220, cal_T_220, cal_T_220, cal_T_220)


    EE_90_90   = (DL_EE.+ssl_EE.+ab_EE.+EE_foregrounds(D_EE_90_90, A_80_EE, α_EE, β_EE,
                                       ν_eff[2,1], ν_eff[2,1], SPT3G_windows_lmax)) ./
                                    _calibration(cal_E_90, cal_E_90, cal_E_90, cal_E_90)

    EE_90_150  = (DL_EE.+ssl_EE.+ab_EE.+EE_foregrounds(D_EE_90_150, A_80_EE, α_EE, β_EE,
                                       ν_eff[2,1], ν_eff[2,2], SPT3G_windows_lmax)) ./
                                    _calibration(cal_E_90, cal_E_150, cal_E_90, cal_E_150)

    EE_90_220  = (DL_EE.+ssl_EE.+ab_EE.+EE_foregrounds(D_EE_90_220, A_80_EE, α_EE, β_EE,
                                       ν_eff[2,1], ν_eff[2,3], SPT3G_windows_lmax)) ./
                                    _calibration(cal_E_90, cal_E_220, cal_E_90, cal_E_220)

    EE_150_150 = (DL_EE.+ssl_EE.+ab_EE.+EE_foregrounds(D_EE_150_150, A_80_EE, α_EE, β_EE,
                                       ν_eff[2,2], ν_eff[2,2], SPT3G_windows_lmax)) ./
                                    _calibration(cal_E_150, cal_E_150, cal_E_150, cal_E_150)

    EE_150_220 = (DL_EE.+ssl_EE.+ab_EE.+EE_foregrounds(D_EE_150_220, A_80_EE, α_EE, β_EE,
                                       ν_eff[2,2], ν_eff[2,3], SPT3G_windows_lmax)) ./
                                    _calibration(cal_E_150, cal_E_220, cal_E_150, cal_E_220)

    EE_220_220 = (DL_EE.+ssl_EE.+ab_EE.+EE_foregrounds(D_EE_220_220, A_80_EE, α_EE, β_EE,
                                       ν_eff[2,3], ν_eff[2,3], SPT3G_windows_lmax)) ./
                                    _calibration(cal_E_220, cal_E_220, cal_E_220, cal_E_220)


    TE_90_90   = (DL_TE.+ssl_TE.+ab_TE.+TE_foregrounds(A_80_TE, α_TE, β_TE,
                                       ν_eff[4,1], ν_eff[4,1], SPT3G_windows_lmax)) ./
                                    _calibration(cal_T_90, cal_E_90, cal_T_90, cal_E_90)

    TE_90_150  = (DL_TE.+ssl_TE.+ab_TE.+TE_foregrounds(A_80_TE, α_TE, β_TE,
                                       ν_eff[4,1], ν_eff[4,2], SPT3G_windows_lmax)) ./
                                    _calibration(cal_T_90, cal_E_150, cal_T_150, cal_E_90)

    TE_90_220  = (DL_TE.+ssl_TE.+ab_TE.+TE_foregrounds(A_80_TE, α_TE, β_TE,
                                       ν_eff[4,1], ν_eff[4,3], SPT3G_windows_lmax)) ./
                                    _calibration(cal_T_90, cal_E_220, cal_T_220, cal_E_90)

    TE_150_150 = (DL_TE.+ssl_TE.+ab_TE.+TE_foregrounds(A_80_TE, α_TE, β_TE,
                                       ν_eff[4,2], ν_eff[4,2], SPT3G_windows_lmax)) ./
                                    _calibration(cal_T_150, cal_E_150, cal_T_150, cal_E_150)

    TE_150_220 = (DL_TE.+ssl_TE.+ab_TE.+TE_foregrounds(A_80_TE, α_TE, β_TE,
                                       ν_eff[4,2], ν_eff[4,3], SPT3G_windows_lmax)) ./
                                    _calibration(cal_T_150, cal_E_220, cal_T_150, cal_E_220)

    TE_220_220 = (DL_TE.+ssl_TE.+ab_TE.+TE_foregrounds(A_80_TE, α_TE, β_TE,
                                       ν_eff[4,3], ν_eff[4,3], SPT3G_windows_lmax)) ./
                                    _calibration(cal_T_220, cal_E_220, cal_T_220, cal_E_220)

    model_matrix = zeros(18, 44)
    model_matrix[1,:]  = window[:,1,:]*TT_90_90
    model_matrix[2,:]  = window[:,2,:]*TE_90_90
    model_matrix[3,:]  = window[:,3,:]*EE_90_90
    model_matrix[4,:]  = window[:,4,:]*TT_90_150
    model_matrix[5,:]  = window[:,5,:]*TE_90_150
    model_matrix[6,:]  = window[:,6,:]*EE_90_150
    model_matrix[7,:]  = window[:,7,:]*TT_90_220
    model_matrix[8,:]  = window[:,8,:]*TE_90_220
    model_matrix[9,:]  = window[:,9,:]*EE_90_220
    model_matrix[10,:] = window[:,10,:]*TT_150_150
    model_matrix[11,:] = window[:,11,:]*TE_150_150
    model_matrix[12,:] = window[:,12,:]*EE_150_150
    model_matrix[13,:] = window[:,13,:]*TT_150_220
    model_matrix[14,:] = window[:,14,:]*TE_150_220
    model_matrix[15,:] = window[:,15,:]*EE_150_220
    model_matrix[16,:] = window[:,16,:]*TT_220_220
    model_matrix[17,:] = window[:,17,:]*TE_220_220
    model_matrix[18,:] = window[:,18,:]*EE_220_220

    return model_matrix
end
