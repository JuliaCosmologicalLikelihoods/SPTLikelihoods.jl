@testset "SPT-3G D1 optional priors" begin
    values = _load_txt_parameters(joinpath(SPT3G_D1_FIXTURE_DIR, "baseline_parameters.txt"))
    parameters = SPT3GD1Parameters(values)
    tau = values["tau"] # 0.054

    # 1. Independent checks of every individual Gaussian prior contribution
    g_pen(val, mean, sigma) = ((val - mean) / sigma)^2 / 2

    p_cib_alpha = g_pen(parameters.cib_alpha, 0.53, 0.1)
    p_dust_ee_amp = g_pen(parameters.dust_ee.amplitude, 0.05, 0.022)
    p_dust_ee_alpha = g_pen(parameters.dust_ee.alpha, -2.42, 0.04)
    p_dust_ee_beta = g_pen(parameters.dust_ee.beta, 1.51, 0.04)
    p_dust_te_amp = g_pen(parameters.dust_te.amplitude, 0.12, 0.051)
    p_dust_te_alpha = g_pen(parameters.dust_te.alpha, -2.42, 0.04)
    p_dust_te_beta = g_pen(parameters.dust_te.beta, 1.51, 0.04)
    p_dust_tt_amp = g_pen(parameters.dust_tt.amplitude, 1.88, 0.96)
    p_dust_tt_alpha = g_pen(parameters.dust_tt.alpha, -2.53, 0.05)
    p_dust_tt_beta = g_pen(parameters.dust_tt.beta, 1.48, 0.02)
    p_tsz = g_pen(parameters.tsz, 3.2279, 2.3764)
    p_ksz = g_pen(parameters.ksz, 3.7287, 4.644)
    p_kappa = g_pen(parameters.kappa, 0.0, 0.00045)
    p_tcal = g_pen(parameters.tcal_ext150, 1.0, 0.00360)
    p_beams = sum(b^2 / 2 for b in parameters.beam_modes)
    p_t2p_1 = g_pen(parameters.t2p[1], -0.006490, 0.001054)
    p_t2p_2 = g_pen(parameters.t2p[2], -0.011941, 0.002103)
    p_t2p_3 = g_pen(parameters.t2p[3], -0.022684, 0.006641)
    p_tau = g_pen(tau, 0.051, 0.006)

    # Verify that the upstream tau prior change (from mean=0.054 to mean=0.051)
    # explains the exact 0.125 discrepancy with the stale upstream test scalar:
    # ((0.054 - 0.051) / 0.006)^2 / 2 = (0.5)^2 / 2 = 0.125.
    @test p_tau ≈ 0.125 atol=1e-15

    expected_nuisance_sum = p_cib_alpha + p_dust_ee_amp + p_dust_ee_alpha + p_dust_ee_beta +
        p_dust_te_amp + p_dust_te_alpha + p_dust_te_beta +
        p_dust_tt_amp + p_dust_tt_alpha + p_dust_tt_beta +
        p_tsz + p_ksz + p_kappa + p_tcal + p_beams +
        p_t2p_1 + p_t2p_2 + p_t2p_3

    @test prior_penalty(parameters) ≈ expected_nuisance_sum rtol=1e-12
    @test prior_penalty(parameters; tau=nothing) ≈ expected_nuisance_sum rtol=1e-12

    expected_total_penalty = expected_nuisance_sum + p_tau
    penalty_with_tau = prior_penalty(parameters; tau=tau)
    @test penalty_with_tau ≈ expected_total_penalty rtol=1e-12
    @test penalty_with_tau ≈ 3.6353251769474317 rtol=1e-12
    @test prior_chi2(parameters; tau=tau) ≈ 2 * penalty_with_tau rtol=1e-14
    @test logprior(parameters; tau=tau) == -penalty_with_tau

    # 2. Reconstruct pinned candl total log-likelihood from data and prior terms
    # candl total = -data_chi2 / 2 - prior_penalty
    data_chi2_baseline = 1373.2802544908752
    reconstructed_candl_total = -data_chi2_baseline / 2 - penalty_with_tau
    @test reconstructed_candl_total ≈ -690.275452422385 rtol=1e-12
    # Verify exact 0.125 shift from the stale upstream scalar -690.15045242
    stale_upstream_loglike = -690.15045242
    @test reconstructed_candl_total ≈ stale_upstream_loglike - 0.125 rtol=1e-8
end
