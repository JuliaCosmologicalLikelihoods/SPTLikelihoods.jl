"""
    test/test_turing.jl

The Turing extension. The likelihood is stated with `~`, not `@addlogprob!`, so
these tests check what that choice buys: the distribution's `logpdf` is this
package's own likelihood plus the Gaussian normalization, `rand` really draws
from the released covariance, the prior table reproduces `prior_penalty`, and
the model's log joint is exactly priors plus likelihood.

This package keeps the Cholesky factor of the covariance, so unlike its sibling
likelihoods it gives up nothing: the density is normalized *and* it can be
sampled from. Both are asserted here, because both are properties of what the
release stores and would be silently lost if that ever changed.
"""

using Turing
using Distributions
using Random

const TURING_EXT = Base.get_extension(SPTLikelihoods, :SPTLikelihoodsTuringExt)

# A small, exactly-known Gaussian, built the way test_spt3g_d1_likelihood.jl
# builds its core case. Defined here rather than reused, so this file does not
# depend on include order.
function synthetic_spt_likelihood()
    ells = collect(2:5)
    windows = [
        [1.0 0.0; 0.0 1.0; 0.0 0.0; 0.0 0.0],
        [0.0 0.0 1.0; 0.0 0.0 0.0; 0.0 1.0 0.0; 1.0 0.0 0.0],
    ]
    data = [1.0, 2.0, 3.0, 4.0, 5.0]
    n = length(data)
    covariance = Matrix(2.0 * I(n))
    for i in 1:n, j in 1:n
        covariance[i, j] += 0.05 * sin(i + j)
    end
    like = SPT3GD1Likelihood(SPT3GD1Data(
        data, covariance, windows, ells;
        spectrum_order=["TT 90x90", "TE 90x90"], bins_per_spectrum=[2, 3],
    ))
    return like, Symmetric(covariance)
end

@testset "SPT-3G D1 Turing extension" begin
    @testset "loaded" begin
        @test TURING_EXT !== nothing
    end

    @testset "the `~` is a normalized Gaussian" begin
        like, covariance = synthetic_spt_likelihood()
        model = [0.9, 1.8, 3.1, 4.2, 4.8]
        d = TURING_EXT.SPT3GD1Bandpowers(like, model)

        @test length(d) == 5

        # Against an explicit Gaussian, not merely a constant offset from one.
        # A dropped, mis-signed, or logdet-of-the-wrong-matrix normalization all
        # fail here. The released 1392 x 1392 case cannot be checked this way.
        reference = MvNormal(model, covariance)
        @test logpdf(d, like.data.data_vector) ≈
            logpdf(reference, like.data.data_vector) rtol=1e-12
        other = [2.0, 1.0, 4.0, 3.0, 6.0]
        @test logpdf(d, other) ≈ logpdf(reference, other) rtol=1e-12

        @test gaussian_normalization(like) ≈
            -5 / 2 * log(2 * pi) - logdet(covariance) / 2 rtol=1e-12
    end

    @testset "rand draws from the released covariance" begin
        # The CamSpec and JointCMB distributions refuse to be sampled, because
        # they store a precision matrix with no factor. This one has the factor,
        # so sampling must work — and must use it correctly.
        like, covariance = synthetic_spt_likelihood()
        model = [0.9, 1.8, 3.1, 4.2, 4.8]
        d = TURING_EXT.SPT3GD1Bandpowers(like, model)

        draw = rand(Random.MersenneTwister(20260920), d)
        @test length(draw) == 5
        @test all(isfinite, draw)

        # Exactly the stored lower factor applied to one normal draw.
        expected = model .+ like.data.cov_chol.L * randn(Random.MersenneTwister(20260920), 5)
        @test draw == expected

        # Independently of that formula: the empirical second moment must be the
        # released covariance, which catches L used where L' was meant.
        rng = Random.MersenneTwister(1234)
        n_draws = 40_000
        draws = reduce(hcat, (rand(rng, d) for _ in 1:n_draws))
        sample_mean = vec(sum(draws; dims=2) ./ n_draws)
        centered = draws .- sample_mean
        empirical = (centered * centered') ./ (n_draws - 1)
        @test maximum(abs, empirical .- covariance) < 0.1
        @test maximum(abs, sample_mean .- model) < 0.05
    end

    @testset "logpdf is the package likelihood" begin
        like = SPT3GD1Likelihood()
        (; params, cmb, model) = _load_test_models()
        prediction = predict(like, model, cmb, params)
        d = TURING_EXT.SPT3GD1Bandpowers(like, prediction)

        @test length(d) == 1392
        # Exactly, not approximately: the distribution runs the same covariance
        # solve `chi2` does, so the only difference from `loglikelihood` is the
        # normalization this package keeps separate.
        @test logpdf(d, like.data.data_vector) ==
            loglikelihood(like, prediction) + gaussian_normalization(like)
        # ... at the pinned baseline chi2 of the multipoint suite.
        @test logpdf(d, like.data.data_vector) ≈
            -1373.2802544908752 / 2 + gaussian_normalization(like) rtol=1e-12

        # A different vector must give a different, finite answer: the
        # distribution must not ignore its argument and echo the data.
        shifted = like.data.data_vector .+ 1.0
        @test isfinite(logpdf(d, shifted))
        @test logpdf(d, shifted) != logpdf(d, like.data.data_vector)

        @test_throws DimensionMismatch TURING_EXT.SPT3GD1Bandpowers(like, prediction[1:10])
    end

    @testset "the prior table restates prior_penalty, not something new" begin
        priors = TURING_EXT.SPT3G_D1_PRIOR_DISTRIBUTIONS
        official = TURING_EXT.SPT3G_D1_GAUSSIAN_PRIORS
        unconstrained = TURING_EXT.SPT3G_D1_UNCONSTRAINED_PARAMETERS

        @test length(priors) == 43
        @test keys(priors) == Symbol.(SPT3G_D1_PARAMETER_NAMES)
        @test length(official) == 26
        @test length(unconstrained) == 17
        # Every parameter is either constrained by the YAML or declared invented.
        # No third category, and no overlap.
        @test isempty(intersect(keys(official), unconstrained))
        @test sort(collect(union(keys(official), unconstrained))) ==
            sort(collect(Symbol.(SPT3G_D1_PARAMETER_NAMES)))

        for name in keys(official)
            @test priors[name] isa Normal
            @test (mean(priors[name]), std(priors[name])) == official[name]
        end

        # The real check. `prior_penalty` is the package's authority; this table
        # is a restatement of it, so the two must differ by the *same* constant
        # at every point. A wrong mean, sigma or missing term breaks that.
        values = _load_txt_parameters(joinpath(SPT3G_D1_FIXTURE_DIR, "baseline_parameters.txt"))
        offsets = Float64[]
        for shift in (0.0, 0.5, -0.3, 1.7)
            point = copy(values)
            for name in keys(official)
                key = string(name)
                point[key] = values[key] + shift * official[name][2]
            end
            table = sum(logpdf(priors[name], point[string(name)]) for name in keys(official))
            push!(offsets, table + prior_penalty(SPT3GD1Parameters(point)))
        end
        @test all(≈(offsets[1], rtol=1e-12), offsets)
        # And that constant is the sum of the omitted log-normalizations.
        expected_offset = -sum(log(official[name][2]) + log(2 * pi) / 2
                               for name in keys(official))
        @test offsets[1] ≈ expected_offset rtol=1e-12
    end

    @testset "log joint is priors plus likelihood" begin
        like = SPT3GD1Likelihood()
        (; params, params_dict, cmb, model) = _load_test_models()
        turing_model = TURING_EXT.spt3g_d1_model(like, model, cmb)

        names = Symbol.(SPT3G_D1_PARAMETER_NAMES)
        values = NamedTuple{names}(
            ntuple(i -> params_dict[SPT3G_D1_PARAMETER_NAMES[i]], length(names)))

        prediction = predict(like, model, cmb, params)
        likelihood_term = loglikelihood(like, prediction) + gaussian_normalization(like)
        priors = TURING_EXT.SPT3G_D1_PRIOR_DISTRIBUTIONS
        prior_term = sum(logpdf(priors[name], values[name]) for name in names)

        @test Turing.DynamicPPL.logjoint(turing_model, values) ≈
            prior_term + likelihood_term rtol=1e-12
        @test isfinite(Turing.DynamicPPL.logjoint(turing_model, values))

        # The log joint matching is itself the ordering check: the 43 `~`
        # statements feed the vector constructor positionally, so any permutation
        # would build different parameters and a different likelihood term.
    end

    @testset "the model is differentiable" begin
        # Establishes that the `~` formulation stays differentiable and that a
        # gradient-based sampler can run. A 2-step NUTS run here does pass, but
        # costs about eight minutes: 43 ForwardDiff duals through the full
        # foreground and instrument model on the 4094-point ell grid. The
        # gradient is the part that would actually break, and test_spt3g_d1_ad.jl
        # already checks its values against Mooncake and finite differences.
        like = SPT3GD1Likelihood()
        (; params_dict, cmb, model) = _load_test_models()
        turing_model = TURING_EXT.spt3g_d1_model(like, model, cmb)

        names = Symbol.(SPT3G_D1_PARAMETER_NAMES)
        start = NamedTuple{names}(
            ntuple(i -> params_dict[SPT3G_D1_PARAMETER_NAMES[i]], length(names)))
        density = Turing.DynamicPPL.LogDensityFunction(
            turing_model, Turing.DynamicPPL.getlogjoint,
            Turing.DynamicPPL.VarInfo(turing_model); adtype=AutoForwardDiff())
        x0 = collect(Float64, start)
        value, gradient = Turing.LogDensityProblems.logdensity_and_gradient(density, x0)
        @test isfinite(value)
        @test length(gradient) == 43
        @test all(isfinite, gradient)
        @test !all(iszero, gradient)
    end
end
