@testset "SPT-3G D1 foreground stages" begin
    fixture = SPT3G_D1_FIXTURE_DIR
    (; params, cmb, model) = _load_test_models(fixture)
    stages = foreground_stages(model, cmb, params)
    fast_pre = SPTLikelihoods._fast_preinstrument_Dls(model, cmb, params)

    # 1. Fast pre-instrument path matches staged foreground pipeline
    @test fast_pre ≈ stages.dust_ee rtol=1e-12 atol=1e-12

    # 2. Cumulative final foreground stage matches candl reference across all 85,974 elements
    ref_preinstrument = vec(_read_matrix(joinpath(fixture, "stage_preinstrument.txt")))
    @test stages.dust_ee ≈ ref_preinstrument rtol=1e-8 atol=1e-8
    @test fast_pre ≈ ref_preinstrument rtol=1e-8 atol=1e-8

    # 3. Individual foreground transformations match candl checkpoints
    chk_stages = (
        :input, :ssl, :aberration, :poisson, :cib, :tsz,
        :tsz_cib, :ksz, :dust_tt, :dust_te, :dust_ee,
    )
    chk = _read_matrix(joinpath(fixture, "foreground_checkpoints.txt"))
    n_ell = length(model.ells)
    for row in 1:size(chk, 1)
        s_idx = Int(chk[row, 1]) + 1
        ell_val = Int(chk[row, 2])
        ell_idx = ell_val - model.ells[1] + 1
        stage_sym = chk_stages[s_idx]
        stage_vec = getproperty(stages, stage_sym)
        for spec in 1:21
            candl_val = chk[row, 2 + spec]
            actual_val = stage_vec[(spec - 1) * n_ell + ell_idx]
            @test actual_val ≈ candl_val rtol=1e-10 atol=1e-10
        end
    end
end
