@testset "SPT-3G D1 candl reference fixture" begin
    reference = read(joinpath(SPT3G_D1_FIXTURE_DIR, "reference.json"), String)
    cmb_mat = _read_matrix(joinpath(SPT3G_D1_FIXTURE_DIR, "cmb_spectra.txt"))
    baseline_unbinned = _read_matrix(joinpath(SPT3G_D1_FIXTURE_DIR, "baseline_final_unbinned.txt"))
    binned_models = _read_matrix(joinpath(SPT3G_D1_FIXTURE_DIR, "binned_models.txt"))
    residuals = _read_matrix(joinpath(SPT3G_D1_FIXTURE_DIR, "residuals.txt"))

    @test occursin("SPT3G_D1_TnE_v0", reference)

    # Provenance must be pinned in the committed fixture, and must match the
    # revisions the loader enforces against the runtime artifact's metadata.json.
    @test occursin(SPTLikelihoods.SPT3G_D1_CANDL_REVISION, reference)
    @test occursin(SPTLikelihoods.SPT3G_D1_SOURCE_REVISION, reference)
    @test occursin("-690.15045242", reference)   # the stale upstream scalar
    # Committed metadata must not record this workstation's checkout location.
    @test !occursin("/home/", reference)
    @test size(cmb_mat) == (4094, 4)
    @test size(baseline_unbinned) == (4094, 21)
    @test size(binned_models) == (1392, 4)
    @test size(residuals) == (1392, 4)
    @test all(isfinite, baseline_unbinned)
    @test all(isfinite, binned_models)
    @test all(isfinite, residuals)

    for name in ("foregrounds", "instrument", "combined")
        mp_unbinned = _read_matrix(joinpath(SPT3G_D1_FIXTURE_DIR, "$(name)_final_unbinned.txt"))
        @test size(mp_unbinned) == (4094, 21)
        @test all(isfinite, mp_unbinned)
    end
end
