@testset "SPT-3G 2018 legacy likelihood smoke test" begin
    ℓs = collect(1:3200)
    DL_TT = fill(1000.0, 3200)
    DL_TE = fill(100.0, 3200)
    DL_EE = fill(50.0, 3200)
    κ = 0.0
    D_TT_90_90 = 10.0; D_TT_90_150 = 8.0; D_TT_90_220 = 15.0
    D_TT_150_150 = 12.0; D_TT_150_220 = 30.0; D_TT_220_220 = 80.0
    D_EE_90_90 = 0.1; D_EE_90_150 = 0.1; D_EE_90_220 = 0.1
    D_EE_150_150 = 0.1; D_EE_150_220 = 0.1; D_EE_220_220 = 0.1
    A_80_cirrus = 2.0; α_cirrus = -2.5; β_cirrus = 1.5
    A_80_cib = 5.0; α_cib = 0.8; β_cib = 1.5
    A_tSZ = 3.0; ξ_tsz_CIB = 0.1; A_kSZ = 2.0
    A_80_EE = 0.05; α_EE = -2.4; β_EE = 1.5
    A_80_TE = 0.10; α_TE = -2.4; β_TE = 1.5
    cal_T_90 = 1.0; cal_T_150 = 1.0; cal_T_220 = 1.0
    cal_E_90 = 1.0; cal_E_150 = 1.0; cal_E_220 = 1.0

    model_matrix = compute_theory(
        DL_TT, DL_TE, DL_EE, κ,
        D_TT_90_90, D_TT_90_150, D_TT_90_220, D_TT_150_150, D_TT_150_220, D_TT_220_220,
        D_EE_90_90, D_EE_90_150, D_EE_90_220, D_EE_150_150, D_EE_150_220, D_EE_220_220,
        A_80_cirrus, α_cirrus, β_cirrus, A_80_cib, α_cib, β_cib, A_tSZ, ξ_tsz_CIB,
        A_kSZ, A_80_EE, α_EE, β_EE, A_80_TE, α_TE, β_TE,
        cal_T_90, cal_T_150, cal_T_220, cal_E_90, cal_E_150, cal_E_220, ℓs,
    )

    @test size(model_matrix) == (18, 44)
    @test all(isfinite, model_matrix)

    res = slice_theory(model_matrix)
    @test length(res) == 728
    @test all(isfinite, res)

    Σ = compute_cov(model_matrix)
    @test size(Σ) == (728, 728)
    @test all(isfinite, Σ)
    @test isposdef(Symmetric(Σ))
end
