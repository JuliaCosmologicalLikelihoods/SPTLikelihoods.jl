# SPTLikelihoods.jl

[![CI](https://github.com/JuliaCosmologicalLikelihoods/SPTLikelihoods.jl/actions/workflows/CI.yml/badge.svg)](https://github.com/JuliaCosmologicalLikelihoods/SPTLikelihoods.jl/actions/workflows/CI.yml)

A Julia-native, fully auto-differentiable library for South Pole Telescope (SPT) CMB likelihoods.

## Supported likelihoods

- **SPT-3G D1 (2019–2020) TT/TE/EE**: standalone, high-performance implementation matching the official `candl` release (43 nuisance parameters, $\ell \in [2, 4095]$, 21 cross-spectra, 1,392 bandpowers).
- **SPT-3G 2018 TT/TE/EE**: pre-existing 2018 likelihood as presented in [Balkenhol et al. (2022)](https://arxiv.org/abs/2212.05642).

Coupled with a differentiable Boltzmann solver or emulator (such as [Capse.jl](https://github.com/CosmologicalEmulators/Capse.jl)), `SPTLikelihoods.jl` supports gradient-based samplers including HMC, NUTS, and Micro-Canonical Hamiltonian Monte Carlo (MCHMC).

---

## Installation

> [!IMPORTANT]
> `SPTLikelihoods.jl` depends on [`CMBForegrounds.jl`](https://github.com/JuliaCosmologicalLikelihoods/CMBForegrounds.jl), which is **not yet registered in the General registry**. Neither package can be installed with `Pkg.add("SPTLikelihoods")` until that happens. Installation by URL is supported today.

**Julia 1.11 and newer** resolve the dependency automatically. `Project.toml` carries a `[sources]` entry pinning `CMBForegrounds` to tag `v0.4.0`, and Pkg follows it recursively:

```julia
using Pkg
Pkg.add(url="https://github.com/JuliaCosmologicalLikelihoods/SPTLikelihoods.jl")
```

**Julia 1.10** does not support `[sources]` and silently ignores it, so the pinned dependency has to be installed first:

```julia
using Pkg
Pkg.add(PackageSpec(
    url="https://github.com/JuliaCosmologicalLikelihoods/CMBForegrounds.jl",
    rev="v0.4.0",
))
Pkg.add(url="https://github.com/JuliaCosmologicalLikelihoods/SPTLikelihoods.jl")
```

No branch is ever tracked: both paths pin the same immutable `v0.4.0` tag. The SPT-3G D1 runtime data is a Julia artifact and is downloaded from Zenodo on first use; see [Data provenance](#data-provenance).

---

## Quickstart: SPT-3G D1

Every Julia block in this section is executed on each CI run by
[`validation/readme_example.jl`](validation/readme_example.jl), so these examples
cannot drift from the code.

```julia
using SPTLikelihoods

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
```

The placeholder spectra above are a smooth stand-in, not a cosmology, so the
chi-square they produce is large by construction. See
[Reproducing the frozen candl reference point](#reproducing-the-frozen-candl-reference-point)
for the pinned release numbers.

Parameters can also be built from a name-keyed `Dict` or `NamedTuple`
(`SPT3GD1Parameters(SPT3G_D1_FIDUCIAL_PARAMETERS)`), and a local converted data
directory can be used instead of the artifact with
`SPT3GD1Likelihood(data_dir)` / `SPT3GD1ForegroundModel(data_dir)`.

### Explicit priors and posterior

`loglikelihood` is **strictly data-only** and never includes implicit priors.
Priors are always evaluated explicitly:

```julia
# Priors are always explicit; `loglikelihood` never includes them.
lp_nuisance = logprior(params)                     # Gaussian nuisance priors
lp_total    = logprior(params; tau=0.054)          # + Planck tau prior
lpost       = logposterior(like, binned_theory, params; tau=0.054)
```

> [!NOTE]
> **$\tau$ prior convention.** The Planck $\tau$ prior is $\tau \sim \mathcal{N}(0.051, 0.006^2)$. At the reference point $\tau = 0.054$ the penalty is
> $\frac{1}{2}\left(\frac{0.054 - 0.051}{0.006}\right)^2 = \frac{1}{2}(0.5)^2 = 0.125$.
> Upstream changed this prior from $0.054 \pm 0.0074$ without updating its stored test scalar, so the published value `-690.15045242` is stale by exactly this `0.125`. This implementation reproduces the corrected total; see below.

### Reproducing the frozen candl reference point

The package ships the frozen `candl` reference spectra and baseline nuisance
parameters as deterministic text fixtures, so the pinned release numbers can be
reproduced directly:

```julia
using SPTLikelihoods
using DelimitedFiles

like  = SPT3GD1Likelihood()
model = SPT3GD1ForegroundModel()

fixture = joinpath(pkgdir(SPTLikelihoods), "test", "fixtures", "spt3g_d1_full")
spectra = readdlm(joinpath(fixture, "cmb_spectra.txt"), Float64; comments=true)
cmb = SPT3GD1CMBTheory(Int.(spectra[:, 1]), spectra[:, 2], spectra[:, 3], spectra[:, 4])

raw = readdlm(joinpath(fixture, "baseline_parameters.txt"), String; comments=true)
params = SPT3GD1Parameters(Dict{String,Float64}(
    raw[i, 1] => parse(Float64, raw[i, 2]) for i in axes(raw, 1)))

binned_theory = predict(like, model, cmb, params)
```

| Quantity | Value |
|---|---:|
| Data $\chi^2$ | `1373.280254490921` |
| Data-only `loglikelihood` | `-686.6401272454605` |
| `logprior(params)` | `-3.510325176947432` |
| `logposterior(...; tau=0.054)` | `-690.2754524224079` |

The test suite pins the same quantities against the frozen `candl` text
fixtures, where the baseline chi-square is `1373.2802544908752`; the artifact
templates reproduce it to a relative `4e-14`. The pinned chi-squares at the
three independent nuisance multipoints, all verified in the test suite, are
`3290.456617577344` (foregrounds), `2134.8981943139897` (instrument) and
`4050.1494799819197` (combined).

The log-posterior above is the corrected total. Upstream's stored scalar
`-690.15045242` differs from it by exactly the `0.125` tau penalty described
in the note above.

---

## Technical specifications (SPT-3G D1)

- **Multipole range**: $\ell \in [2, 4095]$ ($\Delta \ell = 1$, 4,094 multipoles).
- **Frequency channels**: 90 GHz, 150 GHz, 220 GHz.
- **Spectrum ordering (21 spectra)**:
  - 6 TT spectra: `90x90`, `90x150`, `90x220`, `150x150`, `150x220`, `220x220`
  - 9 TE spectra: `90x90`, `90x150`, `150x90`, `90x220`, `220x90`, `150x150`, `150x220`, `220x150`, `220x220` (directional asymmetries preserved)
  - 6 EE spectra: `90x90`, `90x150`, `90x220`, `150x150`, `150x220`, `220x220`
- **Output data vector**: 1,392 bandpowers.
- **Nuisance parameters (43 total)**, in `SPT3G_D1_PARAMETER_NAMES` order:
  - CMB super-sample lensing: `Kappa`
  - Poisson point sources: 6 TT parameters (`TT_Poisson_90x90` … `TT_Poisson_220x220`)
  - CIB clustering: 3 TT parameters (`TT_CIB_150x150`, `TT_CIB_150x220`, `TT_CIB_220x220`) + 1 spectral index (`TT_CIBClustering_Alpha`)
  - Thermal SZ (`TT_tSZ_Amp`), kinetic SZ (`TT_kSZ_Amp`)
  - Galactic cirrus: `TT_GalCirrus_Amp`, `TT_GalCirrus_Alpha`, `TT_GalCirrus_Beta`
  - Polarized Galactic dust: 3 for TE (`TE_PolGalDust_Amp`, `Alpha`, `Beta`) + 3 for EE (`EE_PolGalDust_Amp`, `Alpha`, `Beta`)
  - Temperature-to-polarization leakage: 3 parameters (`T2P2_90`, `T2P2_150`, `T2P2_220`)
  - Main temperature beam eigenmodes: 9 modes (`beta_1` … `beta_9`)
  - Polarization beam transfer: 3 parameters (`beta_pol_90`, `beta_pol_150`, `beta_pol_220`)
  - Calibration: 3 temperature (`Tcal_ext150`, `Tcal_rel90`, `Tcal_rel220`) + 3 polarization (`Ecal_ext150`, `Ecal_rel90`, `Ecal_rel220`)

The 90 GHz channel carries no CIB component by construction. Those legs are
masked structurally before the tSZ–CIB geometric mean, so their CIB derivatives
are exactly zero on both the fast and the staged path. This masking is **not**
extended to a physically active amplitude that merely happens to be zero: the
correlation enters as $\sqrt{A_\mathrm{CIB} D_\mathrm{tSZ}}$, whose slope is
genuinely unbounded as $A_\mathrm{CIB} \to 0$, and the implementation reports
that rather than substituting a convenient zero.

---

## Automatic differentiation

`SPTLikelihoods.jl` provides custom analytical vector-Jacobian products (VJPs) via [ChainRulesCore.jl](https://github.com/JuliaDiff/ChainRulesCore.jl) and supports [Mooncake.jl](https://github.com/compintell/Mooncake.jl) (reverse mode) and [ForwardDiff.jl](https://github.com/JuliaDiff/ForwardDiff.jl) (forward mode) through [DifferentiationInterface.jl](https://github.com/gdalle/DifferentiationInterface.jl).

```julia
using SPTLikelihoods
using DifferentiationInterface
using ADTypes: AutoMooncake
import Mooncake

objective = (x, fixed) -> begin
    like, model, cmb = fixed isa DifferentiationInterface.Constant ? fixed.data : fixed
    params = SPT3GD1Parameters(x)
    return loglikelihood(like, predict(like, model, cmb, params))
end

x0 = Vector(SPT3GD1Parameters())
fixed = DifferentiationInterface.Constant((like, model, cmb))

backend = AutoMooncake(; config=nothing)
prep = prepare_gradient(objective, backend, x0, fixed)
grad = gradient(objective, prep, backend, x0, fixed)
```

Gradients also flow through the **CMB inputs**, which is what makes the package
usable behind a differentiable Boltzmann solver or emulator. The test suite pins
reverse-mode gradients with respect to TT/TE/EE amplitude scalings against
ForwardDiff to better than `1.1e-10`.

---

## Benchmarks

The table below is generated by
[`benchmark/run_benchmarks.jl`](benchmark/run_benchmarks.jl) into
`benchmark/results.md` and is quoted verbatim here and in
[BENCHMARKS.md](BENCHMARKS.md). Hot figures are BenchmarkTools statistics with
compilation excluded; one-shot cold costs are measured once in a fresh process
and are not BenchmarkTools results.

<!-- BEGIN BENCHMARK TABLE -->

<!-- Generated by benchmark/run_benchmarks.jl. Do not edit by hand. -->
<!-- README.md and BENCHMARKS.md must quote this table verbatim. -->

### Environment

| Field | Value |
|---|---|
| Date | 2026-09-17 |
| Julia | 1.12.6 |
| OS | Linux x86_64-linux-gnu, kernel 6.8.0-139-generic |
| CPU | 13th Gen Intel(R) Core(TM) i7-13700H |
| Logical cores | 20 |
| Julia threads | 1 |
| BLAS threads | 10 |
| BenchmarkTools | 1.8.0 |
| CMBForegrounds | 0.4.0 |
| DifferentiationInterface | 0.7.21 |
| ForwardDiff | 1.4.6 |
| Mooncake | 0.5.58 |
| SPTLikelihoods | 0.3.0 |
| Inputs | 21 spectra, 4094 multipoles, 1392 bandpowers, 43 nuisance parameters |

### Hot path (BenchmarkTools, compilation excluded)

| Operation | Median | Minimum | Allocations | Memory |
|---|---:|---:|---:|---:|
| Fast `predict` (production path) | 2.860 ms | 2.597 ms | 1440 | 8.81 MiB |
| Staged `predict` (readable path) | 12.248 ms | 11.531 ms | 1697 | 32.00 MiB |
| Data-only `loglikelihood` | 781.567 μs | 720.857 μs | 6 | 21.89 KiB |
| Combined forward evaluation | 3.955 ms | 3.625 ms | 1448 | 8.84 MiB |
| Hot prepared Mooncake gradient (43 params) | 29.371 ms | 27.966 ms | 165031 | 35.10 MiB |
| Hot prepared ForwardDiff gradient (43 params) | 322.702 ms | 285.243 ms | 6246 | 1.10 GiB |

### One-shot cold costs (fresh process, compilation included)

These are single measurements from `benchmark/cold_prepare.jl`, timed with
`time_ns`/`@timed` and peak RSS growth from `Sys.maxrss`. They are not
BenchmarkTools statistics and have no median or minimum.

| One-shot cost | Wall time | Allocated | Peak RSS growth |
|---|---:|---:|---:|
| Package load (`using SPTLikelihoods`) | 0.366 s | — | 0.0 MiB |
| Artifact-backed construction | 3.021 s | — | — |
| First `predict` (compilation included) | 2.923 s | — | — |
| ForwardDiff `prepare_gradient` | 0.092 s | 12.59 MiB | 13.12 MiB |
| Mooncake `prepare_gradient` (tape build) | 87.355 s | 13106.7 MiB | 551.13 MiB |
| Whole cold process, load to prepared tape | 94.871 s | — | 661.03 MiB |

<!-- END BENCHMARK TABLE -->

Reproduce with:

```bash
julia --project=benchmark -e 'using Pkg; Pkg.instantiate()'
julia --project=benchmark benchmark/run_benchmarks.jl
```

See [BENCHMARKS.md](BENCHMARKS.md) for methodology and the optimisation ledger.

---

## Data provenance

The SPT-3G D1 runtime data is distributed as a Julia artifact, not committed to
this repository. It was converted from the pinned upstream revisions below by
[`validation/convert_spt3g_d1_artifact.py`](validation/convert_spt3g_d1_artifact.py),
and `metadata.json` inside the artifact records the same revision, which the
loader verifies at construction time.

| Field | Value |
|---|---|
| `candl` engine revision | `650db0a6a0a2febed1e57350f0b991b537a4d2f7` |
| `spt_candl_data` revision | `bfe809a140087d19412aa6bc8c8d4ba18840b315` |
| Likelihood descriptor | `SPT3G_D1_TnE_v0/SPT3G_D1_TnE.yaml` |
| Zenodo record | <https://zenodo.org/records/22818351> |
| DOI | `10.5281/zenodo.22818351` |
| Archive | `SPT3G_D1_TnE_v0_20260917.tar.xz` |
| Archive SHA-256 | `2eeabff94fc43c3989bbd57a9328dd3ca097fec3e5123b579fec76dbd587ad6d` |
| Julia artifact tree hash | `e166eed00f35ada69a1f5b6ffb12e853704c49c0` |

---

## Citation

If you use the **SPT-3G D1 T&E** likelihood, cite the references listed by
`spt_candl_data` for this dataset:

- Camphuis et al. (2025), SPT-3G D1 T&E — <https://pole.uchicago.edu/public/Home.html>
- Quan et al. (2025), in preparation
- Balkenhol et al., _`candl`: Cosmic Microwave Background Analysis with a Differentiable Likelihood_ (the `candl` release paper), [arXiv:2401.13433](https://arxiv.org/abs/2401.13433)

If you use the **SPT-3G 2018** likelihood:

- Balkenhol et al. (SPT-3G Collaboration), _Cosmological Constraints from the SPT-3G 2018 TT, TE, and EE Power Spectra_, [arXiv:2212.05642](https://arxiv.org/abs/2212.05642)
- Dutcher et al. (SPT-3G Collaboration), _Measurement of the CMB Temperature and E-Mode Polarization Angular Power Spectra over 1500 Square Degrees of the Southern Sky with SPT-3G_, [arXiv:2101.01684](https://arxiv.org/abs/2101.01684)

If this package is useful in the context of differentiable CMB inference:

- Bonici, Bianchini, and Ruiz-Zapatero, _Capse.jl: efficient and auto-differentiable CMB power spectra emulation_, [arXiv:2307.14339](https://arxiv.org/abs/2307.14339)
