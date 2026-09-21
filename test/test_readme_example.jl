@testset "README example executes" begin
    # `validation/readme_example.jl` is a verbatim copy of the README's Julia
    # blocks. Running it here means a README example that references an undefined
    # name, a renamed constructor, or a removed keyword fails the suite instead of
    # silently rotting in the documentation.
    include(joinpath(@__DIR__, "..", "validation", "readme_example.jl"))

    ctx = READMEExample.quickstart()
    @test length(ctx.binned_theory) == 1392
    @test length(Vector(ctx.params)) == 43
    @test ctx.cmb.ells == collect(2:4095)
    @test all(isfinite, ctx.binned_theory)
    @test isfinite(ctx.chi_squared) && ctx.chi_squared > 0
    @test ctx.log_like ≈ -ctx.chi_squared / 2
    # Priors are strictly additive on top of the data-only log-likelihood
    @test ctx.lp_total ≈ ctx.lp_nuisance - 0.125
    @test ctx.lpost ≈ ctx.log_like + ctx.lp_total

    # The frozen candl reference point reproduces the pinned release numbers.
    ref = READMEExample.reference_point()
    @test ref.chi_squared ≈ 1373.2802544908752 rtol=1e-9
    @test ref.log_like ≈ -ref.chi_squared / 2
    @test ref.lpost ≈ ref.log_like + ref.lp_nuisance - 0.125 rtol=1e-12
    # The stale upstream scalar differs from the correct total by exactly 0.125,
    # which is the tau penalty at tau = 0.054 under N(0.051, 0.006).
    @test ref.lpost ≈ -690.15045242 - 0.125 rtol=1e-8

    ad = READMEExample.ad_example(ctx)
    @test length(ad.x0) == 43
    @test length(ad.grad) == 43
    @test all(isfinite, ad.grad)
end
