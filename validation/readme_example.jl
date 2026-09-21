# Executable copy of the README quickstart.
#
# Every fenced Julia block in README.md under "Quickstart: SPT-3G D1",
# "Explicit priors and posterior", "Reproducing the frozen candl reference
# point", and "Automatic differentiation" appears here verbatim, so a README
# example that stops working fails a test rather than silently rotting.
#
# Run directly:
#     julia --project=. validation/readme_example.jl
#
# It is also executed by `test/test_readme_example.jl` during `Pkg.test()`.

module READMEExample

using SPTLikelihoods
using DelimitedFiles
using DifferentiationInterface
using ADTypes: AutoMooncake
import Mooncake

"""
    quickstart() -> NamedTuple

README section: "Quickstart: SPT-3G D1" and "Explicit priors and posterior".
"""
function quickstart()
    # 1. Artifact-backed data and foreground model. No paths, no downloads to
    #    manage: the SPT3G_D1_TnE_v0 artifact is fetched on first use.
    like  = SPT3GD1Likelihood()
    model = SPT3GD1ForegroundModel()

    # 2. Lensed CMB theory Dℓ in μK², on the D1 grid ℓ = 2:4095. Replace this
    #    placeholder with your Boltzmann solver or emulator output.
    ells = collect(2:4095)
    Dl_TT = @. 5500 * (ells / 220)^(-0.6) * exp(-(ells / 2500)^2) + 60
    Dl_TE = @. 120 * (ells / 500)^0.3 * exp(-(ells / 2000)^2) - 30
    Dl_EE = @. 45 * (ells / 1000)^1.2 * exp(-(ells / 2200)^2) + 0.5
    cmb = SPT3GD1CMBTheory(ells, Dl_TT, Dl_TE, Dl_EE)

    # 3. The 43 nuisance parameters. The zero-argument constructor returns the
    #    fiducial set; SPT3GD1Parameters(x) takes a 43-element vector and
    #    Vector(params) converts back, in SPT3G_D1_PARAMETER_NAMES order.
    params = SPT3GD1Parameters()
    @assert length(Vector(params)) == 43

    # 4. Binned bandpowers (1,392 elements)
    binned_theory = predict(like, model, cmb, params)

    # 5. Data-only chi-square and log-likelihood
    chi_squared = chi2(like, binned_theory)
    log_like    = loglikelihood(like, binned_theory)   # exactly -chi2 / 2

    # 6. Priors are always explicit; `loglikelihood` never includes them.
    lp_nuisance = logprior(params)                     # Gaussian nuisance priors
    lp_total    = logprior(params; tau=0.054)          # + Planck tau prior
    lpost       = logposterior(like, binned_theory, params; tau=0.054)

    return (; like, model, cmb, params, binned_theory,
            chi_squared, log_like, lp_nuisance, lp_total, lpost)
end

"""
    reference_point() -> NamedTuple

README section: "Reproducing the frozen candl reference point". Uses the
deterministic text fixtures shipped with the package.
"""
function reference_point()
    like  = SPT3GD1Likelihood()
    model = SPT3GD1ForegroundModel()

    fixture = joinpath(pkgdir(SPTLikelihoods), "test", "fixtures", "spt3g_d1_full")
    spectra = readdlm(joinpath(fixture, "cmb_spectra.txt"), Float64; comments=true)
    cmb = SPT3GD1CMBTheory(Int.(spectra[:, 1]), spectra[:, 2], spectra[:, 3], spectra[:, 4])

    raw = readdlm(joinpath(fixture, "baseline_parameters.txt"), String; comments=true)
    params = SPT3GD1Parameters(Dict{String,Float64}(
        raw[i, 1] => parse(Float64, raw[i, 2]) for i in axes(raw, 1)))

    binned_theory = predict(like, model, cmb, params)
    return (;
        chi_squared = chi2(like, binned_theory),
        log_like    = loglikelihood(like, binned_theory),
        lp_nuisance = logprior(params),
        lpost       = logposterior(like, binned_theory, params; tau=0.054),
    )
end

"""
    ad_example() -> NamedTuple

README section: "Automatic differentiation".
"""
function ad_example(ctx)
    objective = (x, fixed) -> begin
        like, model, cmb = fixed isa DifferentiationInterface.Constant ? fixed.data : fixed
        params = SPT3GD1Parameters(x)
        return loglikelihood(like, predict(like, model, cmb, params))
    end

    x0 = Vector(SPT3GD1Parameters())
    fixed = DifferentiationInterface.Constant((ctx.like, ctx.model, ctx.cmb))

    backend = AutoMooncake(; config=nothing)
    prep = prepare_gradient(objective, backend, x0, fixed)
    grad = gradient(objective, prep, backend, x0, fixed)

    return (; x0, grad)
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    ctx = READMEExample.quickstart()
    println("binned bandpowers : ", length(ctx.binned_theory))
    println("chi2              : ", ctx.chi_squared)
    println("loglikelihood     : ", ctx.log_like)
    println("logprior          : ", ctx.lp_nuisance)
    println("logprior(tau)     : ", ctx.lp_total)
    println("logposterior      : ", ctx.lpost)

    ref = READMEExample.reference_point()
    println("--- frozen candl reference point ---")
    println("chi2              : ", ref.chi_squared)
    println("loglikelihood     : ", ref.log_like)
    println("logprior          : ", ref.lp_nuisance)
    println("logposterior(tau) : ", ref.lpost)

    ad = READMEExample.ad_example(ctx)
    println("--- prepared Mooncake gradient ---")
    println("length(grad)      : ", length(ad.grad))
    println("all finite        : ", all(isfinite, ad.grad))
end
