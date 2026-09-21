using ChainRulesCore: NoTangent

const _D1_PARAMETER_NAMES = SPT3G_D1_PARAMETER_NAMES

function _d1_parameters_from_vector(x::AbstractVector)
    return SPT3GD1Parameters(x)
end

@testset "SPT-3G D1 custom analytical rrules with arbitrary cotangents" begin
    # 1. _fixed_covariance_solve with deterministic nonuniform arbitrary cotangent
    cov = [2.0 0.5 0.1; 0.5 3.0 0.2; 0.1 0.2 1.5]
    chol = cholesky(Symmetric(cov))
    res = [1.3, -2.1, 0.7]
    ȳ = [0.42, -1.15, 2.31]

    solved, pb_solve = ChainRulesCore.rrule(SPTLikelihoods._fixed_covariance_solve, chol, res)
    dself, dchol, dres = pb_solve(ȳ)

    @test dself == NoTangent()
    @test dchol == NoTangent()

    # Compare residual VJP against ForwardDiff
    fd_dres = ForwardDiff.gradient(r -> dot(ȳ, chol \ r), res)
    @test dres ≈ fd_dres rtol=1e-12

    # 2. _fast_instrument_core with deterministic nonuniform arbitrary cotangent
    fixture = SPT3G_D1_FIXTURE_DIR
    (; model) = _load_test_models(fixture)
    model_s = model
    n_ell_s = length(model_s.ells)

    sky_s = [cos(i * 0.001) for i in 1:(21 * n_ell_s)]
    inst_params_s = [
        -0.007, -0.012, -0.023,
        0.1, -0.2, 0.3, -0.15, 0.05, -0.08, 0.12, -0.04, 0.06,
        0.5, 0.6, 0.7,
        1.002, 0.998, 1.005,
        0.995, 1.003, 0.991,
    ]

    out_s, pb_inst = ChainRulesCore.rrule(SPTLikelihoods._fast_instrument_core, model_s, sky_s, inst_params_s)
    out_cotangent = [sin(i * 0.0005) for i in 1:length(out_s)]
    dself_i, dmodel_i, dsky_s, dparams_s = pb_inst(out_cotangent)

    @test dself_i == NoTangent()
    @test dmodel_i == NoTangent()

    # Test all 21 instrument parameters VJP against ForwardDiff
    fd_dparams = ForwardDiff.gradient(
        p -> dot(out_cotangent, SPTLikelihoods._fast_instrument_core(model_s, sky_s, p)),
        inst_params_s,
    )
    @test dparams_s ≈ fd_dparams rtol=1e-10 atol=1e-10

    # Test sky VJP along deterministic non-uniform direction vectors
    v_sky = [sin(i * 0.01) for i in 1:(21 * n_ell_s)]
    dir_deriv_sky = ForwardDiff.derivative(
        t -> dot(out_cotangent, SPTLikelihoods._fast_instrument_core(model_s, sky_s .+ t .* v_sky, inst_params_s)),
        0.0,
    )
    @test dot(dsky_s, v_sky) ≈ dir_deriv_sky rtol=1e-10 atol=1e-10
end

@testset "SPT-3G D1 inactive 90-GHz cross-spectra and CIB derivative semantics" begin
    fixture = SPT3G_D1_FIXTURE_DIR
    (; params, params_dict, cmb, model) = _load_test_models(fixture)
    values = params_dict
    x0 = [values[name] for name in _D1_PARAMETER_NAMES]

    # Verify that inactive 90-GHz leg is handled structurally before sqrt
    # In SPT3G_D1_SPECTRUM_ORDER: TT 90x90 is index 1, TT 90x150 is index 4, TT 90x220 is index 8
    n_ell = length(model.ells)

    obj_90x90 = function (x, fixed)
        m, c = fixed
        p = SPT3GD1Parameters(x)
        pre = SPTLikelihoods._fast_preinstrument_Dls(m, c, p)
        return sum(@view(pre[1:n_ell]))
    end

    obj_90x150 = function (x, fixed)
        m, c = fixed
        p = SPT3GD1Parameters(x)
        pre = SPTLikelihoods._fast_preinstrument_Dls(m, c, p)
        return sum(@view(pre[(4 - 1) * n_ell + 1:4 * n_ell]))
    end

    obj_90x220 = function (x, fixed)
        m, c = fixed
        p = SPT3GD1Parameters(x)
        pre = SPTLikelihoods._fast_preinstrument_Dls(m, c, p)
        return sum(@view(pre[(8 - 1) * n_ell + 1:8 * n_ell]))
    end

    fixed = DifferentiationInterface.Constant((model, cmb))

    for (name, obj) in (("90x90", obj_90x90), ("90x150", obj_90x150), ("90x220", obj_90x220))
        g_fd = DifferentiationInterface.gradient(obj, AutoForwardDiff(), x0, fixed)
        @test all(isfinite, g_fd)
    end

    # Verify Mooncake gradient matches ForwardDiff on cross-spectrum
    prep_90x150 = DifferentiationInterface.prepare_gradient(obj_90x150, AutoMooncake(; config=nothing), x0, fixed)
    g_mc_90x150 = DifferentiationInterface.gradient(obj_90x150, prep_90x150, AutoMooncake(; config=nothing), x0, fixed)
    g_fd_90x150 = DifferentiationInterface.gradient(obj_90x150, AutoForwardDiff(), x0, fixed)
    @test all(isfinite, g_mc_90x150)
    @test g_mc_90x150 ≈ g_fd_90x150 rtol=1e-8 atol=1e-8

    # Specific semantic checks:
    # 1. TT 90x90 has NO CIB or CIB-tSZ correlation, so gradient w.r.t. all CIB amplitudes must be zero
    g_90x90 = DifferentiationInterface.gradient(obj_90x90, AutoForwardDiff(), x0, fixed)
    # Indices 8, 9, 10 correspond to TT_CIB_150x150, TT_CIB_150x220, TT_CIB_220x220
    @test g_90x90[8] == 0.0
    @test g_90x90[9] == 0.0
    @test g_90x90[10] == 0.0

    # 2. TT 90x150 has active 150 GHz CIB leg, so gradient w.r.t. TT_CIB_150x150 is non-zero
    @test g_fd_90x150[8] != 0.0
    @test isfinite(g_fd_90x150[8])
    # But has NO 220 GHz leg, so TT_CIB_220x220 is zero
    @test g_fd_90x150[10] == 0.0

    # 3. TT 90x220 has active 220 GHz CIB leg, so gradient w.r.t. TT_CIB_220x220 is non-zero
    g_90x220 = DifferentiationInterface.gradient(obj_90x220, AutoForwardDiff(), x0, fixed)
    @test g_90x220[10] != 0.0
    @test isfinite(g_90x220[10])
    # But has NO 150 GHz leg, so TT_CIB_150x150 is zero
    @test g_90x220[8] == 0.0
end

@testset "SPT-3G D1 complete public likelihood differentiation" begin
    fixture = SPT3G_D1_FIXTURE_DIR
    (; params, params_dict, cmb, model) = _load_test_models(fixture)
    values = params_dict
    x0 = [values[name] for name in _D1_PARAMETER_NAMES]

    # Full likelihood loaded directly from bound official artifact (runs unconditionally)
    like = SPT3GD1Likelihood()

    # 1. Differentiate complete public data log-likelihood through windows and covariance
    full_objective = function (x, fixed)
        ctx = fixed isa DifferentiationInterface.Constant ? fixed.data : fixed
        l, m, c = ctx
        p = SPT3GD1Parameters(x)
        pred = predict(l, m, c, p)
        return loglikelihood(l, pred)
    end
    fixed = DifferentiationInterface.Constant((like, model, cmb))

    forwarddiff = AutoForwardDiff()
    mooncake = AutoMooncake(; config=nothing)

    g_fd = DifferentiationInterface.gradient(full_objective, forwarddiff, x0, fixed)
    @test all(isfinite, g_fd)

    prep = DifferentiationInterface.prepare_gradient(full_objective, mooncake, x0, fixed)
    g_mc = DifferentiationInterface.gradient(full_objective, prep, mooncake, x0, fixed)
    @test all(isfinite, g_mc)
    @test g_mc ≈ g_fd rtol=1e-5 atol=1e-5

    # 2. Preparation reuse with changed parameters (guards against constant capture)
    x1 = copy(x0)
    x1[1] += 1e-5     # change Kappa
    x1[2] += 0.5      # change Poisson
    x1[23] += 1e-4    # change T2P

    val0 = full_objective(x0, fixed)
    val1 = full_objective(x1, fixed)
    @test val0 != val1

    g_mc1 = DifferentiationInterface.gradient(full_objective, prep, mooncake, x1, fixed)
    g_fd1 = DifferentiationInterface.gradient(full_objective, forwarddiff, x1, fixed)
    @test all(isfinite, g_mc1)
    @test g_mc1 ≈ g_fd1 rtol=1e-5 atol=1e-5
    @test g_mc1 != g_mc

    # 3. Finite differences directional check along a nontrivial parameter direction
    finite = AutoFiniteDifferences(; fdm=central_fdm(5, 1))
    v_dir = normalize(ones(length(x0)))
    directional_obj = t -> full_objective(x0 .+ t .* v_dir, fixed)
    dir_deriv_fd = dot(g_fd, v_dir)
    dir_deriv_finite = DifferentiationInterface.derivative(directional_obj, finite, 0.0)
    @test isapprox(dir_deriv_finite, dir_deriv_fd, rtol=1e-4)

    # 4. Differentiate logprior separately
    prior_obj = x -> logprior(SPT3GD1Parameters(x))
    g_prior_fd = DifferentiationInterface.gradient(prior_obj, forwarddiff, x0)
    g_prior_mc = DifferentiationInterface.gradient(prior_obj, mooncake, x0)
    @test all(isfinite, g_prior_fd)
    @test g_prior_mc ≈ g_prior_fd rtol=1e-6 atol=1e-7

    # 5. Differentiate full posterior gradient (with tau)
    tau_val = values["tau"]
    posterior_obj = function (x, fixed)
        ctx = fixed isa DifferentiationInterface.Constant ? fixed.data : fixed
        l, m, c = ctx
        p = SPT3GD1Parameters(x)
        pred = predict(l, m, c, p)
        return logposterior(l, pred, p; tau=tau_val)
    end
    g_post_fd = DifferentiationInterface.gradient(posterior_obj, forwarddiff, x0, fixed)
    prep_post = DifferentiationInterface.prepare_gradient(posterior_obj, mooncake, x0, fixed)
    g_post_mc = DifferentiationInterface.gradient(posterior_obj, prep_post, mooncake, x0, fixed)
    @test all(isfinite, g_post_fd)
    @test g_post_mc ≈ g_post_fd rtol=1e-5 atol=1e-5
    @test g_post_fd ≈ g_fd .+ g_prior_fd rtol=1e-10
end

@testset "SPT-3G D1 staged foreground differentiation" begin
    fixture = SPT3G_D1_FIXTURE_DIR
    (; params, params_dict, cmb, model) = _load_test_models(fixture)
    x0 = SPTLikelihoods._to_vector(params)
    n_ell = length(model.ells)
    mooncake = AutoMooncake(; config=nothing)

    # The exported staged path must be differentiable, not only the fast
    # production path. Regression guard for the staged CIB-tSZ NaN defect, where
    # `_tsz_cib_long` took square roots through structurally inactive 90-GHz CIB
    # legs and poisoned all 43 entries of every staged gradient.
    for stage_sym in (:tsz_cib, :dust_tt, :dust_te, :dust_ee)
        f_stage = x -> sum(getfield(foreground_stages(model, cmb, SPT3GD1Parameters(x)), stage_sym))
        g_fd = ForwardDiff.gradient(f_stage, x0)
        @test length(g_fd) == 43
        @test count(isnan, g_fd) == 0
        @test all(isfinite, g_fd)
    end

    # 2. Compare staged dust_ee gradient against fast pre-instrument gradient
    f_staged = x -> sum(foreground_stages(model, cmb, SPT3GD1Parameters(x)).dust_ee)
    f_fast   = x -> sum(SPTLikelihoods._fast_preinstrument_Dls(model, cmb, SPT3GD1Parameters(x)))

    g_staged_fd = ForwardDiff.gradient(f_staged, x0)
    g_fast_fd   = ForwardDiff.gradient(f_fast, x0)
    @test all(isfinite, g_staged_fd)
    @test g_staged_fd ≈ g_fast_fd rtol=1e-10 atol=1e-10

    # 3. Prepared Mooncake agreement on staged foreground output
    staged_obj = (x, fixed) -> begin
        m, c = fixed isa DifferentiationInterface.Constant ? fixed.data : fixed
        sum(foreground_stages(m, c, SPT3GD1Parameters(x)).dust_ee)
    end
    fixed_mc = DifferentiationInterface.Constant((model, cmb))
    prep_staged = DifferentiationInterface.prepare_gradient(staged_obj, mooncake, x0, fixed_mc)
    g_staged_mc = DifferentiationInterface.gradient(staged_obj, prep_staged, mooncake, x0, fixed_mc)
    @test all(isfinite, g_staged_mc)
    @test g_staged_mc ≈ g_staged_fd rtol=1e-6 atol=1e-6

    # 4. Inactive 90-GHz CIB legs on the STAGED path specifically. The fast path
    #    is covered separately; this is the path that previously produced NaN.
    #    `_D1_TT_PAIRS` fixes leg 1 (90 GHz) as carrying no CIB, so CIB gradients
    #    of any 90-GHz auto spectrum are structurally zero, while a cross with an
    #    active leg is finite and non-zero.
    cib_indices = (
        findfirst(==("TT_CIB_150x150"), _D1_PARAMETER_NAMES),
        findfirst(==("TT_CIB_150x220"), _D1_PARAMETER_NAMES),
        findfirst(==("TT_CIB_220x220"), _D1_PARAMETER_NAMES),
    )
    staged_block = (x, spectrum_index) -> sum(@view(
        foreground_stages(model, cmb, SPT3GD1Parameters(x)).tsz_cib[
            (spectrum_index - 1) * n_ell + 1 : spectrum_index * n_ell]))

    g_stage_90x90 = ForwardDiff.gradient(x -> staged_block(x, 1), x0)
    @test all(isfinite, g_stage_90x90)
    for index in cib_indices
        @test g_stage_90x90[index] == 0.0
    end

    g_stage_90x150 = ForwardDiff.gradient(x -> staged_block(x, 4), x0)
    @test all(isfinite, g_stage_90x150)
    @test g_stage_90x150[cib_indices[1]] != 0.0
    @test g_stage_90x150[cib_indices[3]] == 0.0

    g_stage_90x220 = ForwardDiff.gradient(x -> staged_block(x, 8), x0)
    @test all(isfinite, g_stage_90x220)
    @test g_stage_90x220[cib_indices[3]] != 0.0
    @test g_stage_90x220[cib_indices[1]] == 0.0

    # 5. Zero-amplitude boundary policy for PHYSICALLY ACTIVE legs.
    #
    #    Policy: structural masking applies only to legs that carry no CIB by
    #    construction (every 90-GHz leg). It is NOT extended to an active leg
    #    whose amplitude merely happens to evaluate to zero. The tSZ-CIB
    #    correlation enters as sqrt(A_CIB * D_tSZ), whose slope in A_CIB is
    #    genuinely unbounded as A_CIB -> 0, so the honest derivative there is
    #    non-finite. The implementation reports that rather than substituting a
    #    convenient zero, and the staged and fast paths agree on it exactly.
    x_zero_cib = copy(x0)
    x_zero_cib[cib_indices[1]] = 0.0
    tsz_index = findfirst(==("TT_tSZ_Amp"), _D1_PARAMETER_NAMES)

    fast_block_150 = x -> sum(@view(
        SPTLikelihoods._fast_preinstrument_Dls(model, cmb, SPT3GD1Parameters(x))[
            (12 - 1) * n_ell + 1 : 12 * n_ell]))
    g_zero_staged = ForwardDiff.gradient(x -> staged_block(x, 12), x_zero_cib)
    g_zero_fast   = ForwardDiff.gradient(fast_block_150, x_zero_cib)

    @test !isfinite(g_zero_staged[cib_indices[1]])
    @test g_zero_staged[cib_indices[1]] == -Inf
    @test isnan(g_zero_staged[tsz_index])
    # Staged and fast paths must agree on the boundary, not merely away from it
    @test g_zero_fast[cib_indices[1]] == g_zero_staged[cib_indices[1]]
    @test isnan(g_zero_fast[tsz_index]) == isnan(g_zero_staged[tsz_index])

    # Away from the boundary the same active amplitude is ordinary and finite
    g_active = ForwardDiff.gradient(x -> staged_block(x, 12), x0)
    @test all(isfinite, g_active)
    @test g_active[cib_indices[1]] != 0.0
end

@testset "SPT-3G D1 reverse-mode differentiation through CMB inputs" begin
    fixture = SPT3G_D1_FIXTURE_DIR
    (; params, cmb, model) = _load_test_models(fixture)
    like = SPT3GD1Likelihood()

    # Differentiate data log-likelihood w.r.t. CMB amplitude scaling [A_TT, A_TE, A_EE]
    cmb_scale_obj = function (a, fixed)
        ctx = fixed isa DifferentiationInterface.Constant ? fixed.data : fixed
        l, m, c, p = ctx
        scaled_cmb = SPT3GD1CMBTheory(
            c.ells,
            a[1] .* c.TT,
            a[2] .* c.TE,
            a[3] .* c.EE,
        )
        pred = predict(l, m, scaled_cmb, p)
        return loglikelihood(l, pred)
    end
    fixed = DifferentiationInterface.Constant((like, model, cmb, params))
    a0 = [1.0, 1.0, 1.0]

    forwarddiff = AutoForwardDiff()
    mooncake = AutoMooncake(; config=nothing)

    g_fd = DifferentiationInterface.gradient(cmb_scale_obj, forwarddiff, a0, fixed)
    @test all(isfinite, g_fd)

    prep = DifferentiationInterface.prepare_gradient(cmb_scale_obj, mooncake, a0, fixed)
    g_mc = DifferentiationInterface.gradient(cmb_scale_obj, prep, mooncake, a0, fixed)
    @test all(isfinite, g_mc)
    @test g_mc ≈ g_fd rtol=1e-8 atol=1e-8
    @test maximum(abs, g_mc .- g_fd) < 1e-8

    # Pinned reference gradient. Differentiating through the CMB inputs is the
    # advertised coupling to differentiable Boltzmann solvers and emulators, so
    # the values themselves are frozen here, not only ForwardDiff/Mooncake
    # agreement: a change that moved both backends together would otherwise pass.
    g_reference = [-0.6253048029325117, -110.89179744330002, 145.72695961808924]
    @test g_fd ≈ g_reference rtol=1e-10 atol=1e-10
    @test g_mc ≈ g_reference rtol=1e-8 atol=1e-8

    # Reusing the preparation must track a changed input rather than capture it
    a1 = [1.05, 0.97, 1.02]
    g_mc_moved = DifferentiationInterface.gradient(cmb_scale_obj, prep, mooncake, a1, fixed)
    g_fd_moved = DifferentiationInterface.gradient(cmb_scale_obj, forwarddiff, a1, fixed)
    @test all(isfinite, g_mc_moved)
    @test !(g_mc_moved ≈ g_mc)
    @test g_mc_moved ≈ g_fd_moved rtol=1e-8 atol=1e-8
end

