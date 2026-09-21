@testset "SPT-3G D1 Gaussian core" begin
    ells = collect(2:5)
    windows = [
        [1.0 0.0; 0.0 1.0; 0.0 0.0; 0.0 0.0],
        [0.0 0.0 1.0; 0.0 0.0 0.0; 0.0 1.0 0.0; 1.0 0.0 0.0],
    ]
    data = [1.0, 2.0, 3.0, 4.0, 5.0]
    covariance = Matrix{Float64}(I, length(data), length(data))
    fixed_data = SPT3GD1Data(
        data,
        covariance,
        windows,
        ells;
        spectrum_order=["TT 90x90", "TE 90x90"],
        bins_per_spectrum=[2, 3],
    )
    like = SPT3GD1Likelihood(fixed_data)
    unbinned = [1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0]
    model = bin_theory(like, unbinned)

    @test model == [1.0, 2.0, 8.0, 7.0, 5.0]
    @test chi2(like, model) == 34.0
    @test chi2(like, unbinned, Val(:unbinned)) == 34.0
    @test loglikelihood(like, model) == -17.0
    @test loglikelihood(like, unbinned, Val(:unbinned)) == -17.0
    @test_throws DimensionMismatch bin_theory(like, ones(7))
    @test_throws DimensionMismatch chi2(like, ones(4))

    # Test grid mismatch validation in predict
    cmb_mismatch = SPT3GD1CMBTheory(collect(3:6), ones(4), ones(4), ones(4))
    fg_model_match = SPT3GD1ForegroundModel(
        collect(2:4000), ones(3999), ones(3999);
        beam_modes=zeros(3, 3999, 9), main_temperature_beams=ones(3, 3999),
    )
    mock_params = SPT3GD1Parameters(Dict(
        "Kappa" => 0.0, "TT_Poisson_90x90" => 0.0, "TT_Poisson_90x150" => 0.0,
        "TT_Poisson_90x220" => 0.0, "TT_Poisson_150x150" => 0.0, "TT_Poisson_150x220" => 0.0,
        "TT_Poisson_220x220" => 0.0, "TT_CIB_150x150" => 0.0, "TT_CIB_150x220" => 0.0,
        "TT_CIB_220x220" => 0.0, "TT_CIBClustering_Alpha" => 0.0, "TT_tSZ_Amp" => 0.0,
        "TT_kSZ_Amp" => 0.0, "TT_GalCirrus_Amp" => 0.0, "TT_GalCirrus_Alpha" => 0.0,
        "TT_GalCirrus_Beta" => 0.0, "TE_PolGalDust_Amp" => 0.0, "TE_PolGalDust_Alpha" => 0.0,
        "TE_PolGalDust_Beta" => 0.0, "EE_PolGalDust_Amp" => 0.0, "EE_PolGalDust_Alpha" => 0.0,
        "EE_PolGalDust_Beta" => 0.0, "T2P2_90" => 0.0, "T2P2_150" => 0.0, "T2P2_220" => 0.0,
        "beta_1" => 0.0, "beta_2" => 0.0, "beta_3" => 0.0, "beta_4" => 0.0, "beta_5" => 0.0,
        "beta_6" => 0.0, "beta_7" => 0.0, "beta_8" => 0.0, "beta_9" => 0.0,
        "beta_pol_90" => 0.0, "beta_pol_150" => 0.0, "beta_pol_220" => 0.0,
        "Tcal_ext150" => 1.0, "Tcal_rel90" => 1.0, "Tcal_rel220" => 1.0,
        "Ecal_ext150" => 1.0, "Ecal_rel90" => 1.0, "Ecal_rel220" => 1.0,
    ))
    @test_throws ArgumentError predict(like, fg_model_match, cmb_mismatch, mock_params)

    # Core function inference tests
    @inferred chi2(like, model)
    @inferred loglikelihood(like, model)
    @inferred bin_theory(like, unbinned)
    @inferred SPTLikelihoods._fixed_covariance_solve(like.data.cov_chol, [1.0, 2.0, 3.0, 4.0, 5.0])

    # Zero-argument constructors load the bound artifact
    auto_like = SPT3GD1Likelihood()
    @test length(auto_like.data.data_vector) == 1392
    @test length(auto_like.data.windows) == 21
    @inferred chi2(auto_like, auto_like.data.data_vector)
    @inferred loglikelihood(auto_like, auto_like.data.data_vector)

    auto_model = SPT3GD1ForegroundModel()
    @test auto_model.ells == collect(2:4095)

    # Validate rejection of corrupt or mismatched metadata.json
    mktempdir() do tmp
        write(joinpath(tmp, "data_vector.txt"), "1.0\n2.0\n")
        write(joinpath(tmp, "metadata.json"), "{\"spectrum_order\": [\"CORRUPTED_ORDER\"]}")
        @test_throws ArgumentError SPTLikelihoods.load_spt3g_d1_data(tmp)
    end
end

