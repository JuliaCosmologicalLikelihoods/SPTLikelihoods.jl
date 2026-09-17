@testset "SPT-3G D1 instrument stages" begin
    fixture = SPT3G_D1_FIXTURE_DIR
    (; params, cmb, model) = _load_test_models(fixture)
    foreground = foreground_stages(model, cmb, params)
    stages = instrument_stages(model, foreground.dust_ee, params)
    fast_preinstrument = SPTLikelihoods._fast_preinstrument_Dls(model, cmb, params)
    fast_calibration = SPTLikelihoods._fast_instrument_Dls(model, fast_preinstrument, params)

    ref_pre = vec(_read_matrix(joinpath(fixture, "stage_preinstrument.txt")))
    ref_leak = vec(_read_matrix(joinpath(fixture, "stage_leakage.txt")))
    ref_beam = vec(_read_matrix(joinpath(fixture, "stage_beam.txt")))
    ref_cal = vec(_read_matrix(joinpath(fixture, "baseline_final_unbinned.txt")))

    @test fast_preinstrument ≈ ref_pre rtol=1e-8 atol=1e-8
    @test stages.leakage ≈ ref_leak rtol=1e-8 atol=1e-8
    @test stages.beam ≈ ref_beam rtol=1e-8 atol=1e-8
    @test stages.calibration ≈ ref_cal rtol=1e-8 atol=1e-8
    @test fast_calibration ≈ stages.calibration rtol=1e-12 atol=1e-12
end

