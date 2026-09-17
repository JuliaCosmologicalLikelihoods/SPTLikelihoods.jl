@testset "SPT-3G D1 multipoint forward parity and public predict" begin
    fixture = SPT3G_D1_FIXTURE_DIR
    (; params, params_dict, cmb, model) = _load_test_models(fixture)
    multipoint_values = _load_txt_parameters(joinpath(fixture, "multipoint_parameters.txt"))
    baseline = copy(params_dict)

    pinned_chi2 = Dict(
        "baseline" => 1373.2802544908752,
        "foregrounds" => 3290.456617577344,
        "instrument" => 2134.8981943139897,
        "combined" => 4050.1494799819197,
    )

    case_indices = Dict("baseline" => 1, "foregrounds" => 2, "instrument" => 3, "combined" => 4)
    ref_binned_models = _read_matrix(joinpath(fixture, "binned_models.txt"))
    ref_residuals = _read_matrix(joinpath(fixture, "residuals.txt"))

    # Test likelihood from test fixture (always available, runs unconditionally)
    test_like = SPT3GD1Likelihood(SPT3G_D1_TEST_LIKE_DIR)

    # Full likelihood loaded directly from bound official artifact (runs unconditionally)
    real_like = SPT3GD1Likelihood()

    for name in ("baseline", "foregrounds", "instrument", "combined")
        values = copy(baseline)
        if name != "baseline"
            prefix = "$(name)__"
            for (key, value) in multipoint_values
                startswith(key, prefix) || continue
                values[key[length(prefix)+1:end]] = value
            end
        end
        parameters = SPT3GD1Parameters(values)

        # 1. Staged path
        foreground = foreground_stages(model, cmb, parameters)
        instrument = instrument_stages(model, foreground.dust_ee, parameters)
        unbinned_staged = instrument.calibration

        # 2. Fast path
        fast_pre = SPTLikelihoods._fast_preinstrument_Dls(model, cmb, parameters)
        fast_unbinned = SPTLikelihoods._fast_instrument_Dls(model, fast_pre, parameters)

        # Staged and fast paths must agree for all multipoint cases
        @test fast_unbinned ≈ unbinned_staged rtol=1e-12 atol=1e-12

        # Final unbinned vector must agree with frozen candl reference
        ref_unbinned = vec(_read_matrix(joinpath(fixture, "$(name)_final_unbinned.txt")))
        @test fast_unbinned ≈ ref_unbinned rtol=1e-8 atol=1e-8

        # 3. Public likelihood evaluation on test fixture (unconditional)
        test_pred = predict(test_like, model, cmb, parameters)
        @test test_pred ≈ bin_theory(test_like, fast_unbinned) rtol=1e-12 atol=1e-12
        test_c2 = chi2(test_like, test_pred)
        @test isfinite(test_c2)
        test_ll = loglikelihood(test_like, test_pred)
        @test test_ll == -test_c2 / 2

        # 4. Full public path tests with real data
        pred = predict(real_like, model, cmb, parameters)
        @test pred ≈ bin_theory(real_like, fast_unbinned) rtol=1e-12 atol=1e-12

        col = case_indices[name]
        ref_binned = ref_binned_models[:, col]
        @test pred ≈ ref_binned rtol=1e-10 atol=1e-10

        residual = real_like.data.data_vector .- pred
        ref_res = ref_residuals[:, col]
        @test residual ≈ ref_res rtol=1e-10 atol=1e-10

        c2 = chi2(real_like, pred)
        @test c2 ≈ pinned_chi2[name] rtol=1e-9 atol=1e-8

        ll = loglikelihood(real_like, pred)
        @test ll ≈ -pinned_chi2[name] / 2 rtol=1e-9 atol=1e-8
        @test ll == -c2 / 2
    end

    # 5. Verify directional TE spectra remain distinct through leakage, beams, and calibration
    inst_values = copy(baseline)
    for (key, value) in multipoint_values
        startswith(key, "instrument__") || continue
        inst_values[key[length("instrument__")+1:end]] = value
    end
    inst_params = SPT3GD1Parameters(inst_values)
    inst_pre = SPTLikelihoods._fast_preinstrument_Dls(model, cmb, inst_params)
    inst_unbinned = SPTLikelihoods._fast_instrument_Dls(model, inst_pre, inst_params)
    n_ell = length(model.ells)
    te_90_150 = inst_unbinned[(5 - 1) * n_ell + 1:5 * n_ell]
    te_150_90 = inst_unbinned[(6 - 1) * n_ell + 1:6 * n_ell]
    @test maximum(abs.(te_90_150 .- te_150_90)) > 1.0

    binned_inst = predict(real_like, model, cmb, inst_params)
    b_90_150 = binned_inst[249:320]
    b_150_90 = binned_inst[321:392]
    @test maximum(abs.(b_90_150 .- b_150_90)) > 0.1
end
