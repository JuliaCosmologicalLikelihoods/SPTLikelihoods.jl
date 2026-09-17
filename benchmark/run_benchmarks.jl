# Reproducible SPT-3G D1 benchmark harness.
#
#     julia --project=benchmark -e 'using Pkg; Pkg.instantiate()'
#     julia --project=benchmark benchmark/run_benchmarks.jl
#
# Two kinds of measurement, kept strictly apart:
#
#   * Hot, repeatable operations are measured with BenchmarkTools `@benchmark`,
#     with `$` interpolation so setup is never timed. Stateful prepared-gradient
#     calls use `evals=1` so a mutated cache is never reused inside one sample.
#   * One-shot cold costs (package load, artifact construction, first compiled
#     call, AD preparation, peak RSS growth) cannot be sampled repeatedly in a
#     warm process, so they are measured once in a fresh child process by
#     `benchmark/cold_prepare.jl`. They are reported in their own table and are
#     never described as BenchmarkTools results.
#
# Writes the single authoritative results table to `benchmark/results.md`.
# README.md and BENCHMARKS.md quote that file; they are not edited by hand.

using SPTLikelihoods
using LinearAlgebra
using DelimitedFiles
using BenchmarkTools
using ADTypes: AutoForwardDiff, AutoMooncake
using DifferentiationInterface
using Pkg
import Mooncake
import ForwardDiff

const RESULTS_PATH = joinpath(@__DIR__, "results.md")

_kernel_release() = try
    Sys.isunix() ? readchomp(`uname -r`) : "unknown"
catch
    "unknown"
end

function environment_report()
    cpu = Sys.cpu_info()
    deps = Pkg.dependencies()
    versions = Dict{String,String}()
    for (_, info) in deps
        if info.name in ("BenchmarkTools", "DifferentiationInterface", "Mooncake",
                         "ForwardDiff", "CMBForegrounds", "SPTLikelihoods")
            versions[info.name] = string(info.version)
        end
    end
    return (
        julia = string(VERSION),
        os = "$(Sys.KERNEL) $(Sys.MACHINE), kernel $(_kernel_release())",
        cpu = cpu[1].model,
        cores = length(cpu),
        threads = Threads.nthreads(),
        blas_threads = BLAS.get_num_threads(),
        versions = versions,
        date = Base.Libc.strftime("%Y-%m-%d", time()),
    )
end

function run_cold_harness()
    script = joinpath(@__DIR__, "cold_prepare.jl")
    project = Base.active_project()
    output = read(`$(Base.julia_cmd()) --project=$project $script`, String)
    cold = Dict{String,Float64}()
    for line in eachsplit(strip(output), '\n')
        parts = split(line, '=')
        length(parts) == 2 || continue
        cold[String(parts[1])] = parse(Float64, parts[2])
    end
    return cold
end

# ---------------------------------------------------------------- setup ------
env = environment_report()
println("=== SPT-3G D1 benchmark environment ===")
println("Julia   : ", env.julia)
println("OS      : ", env.os)
println("CPU     : ", env.cpu, " ($(env.cores) logical cores)")
println("Threads : ", env.threads, " Julia, ", env.blas_threads, " BLAS")
for (name, version) in sort(collect(env.versions))
    println("  ", rpad(name, 26), version)
end

fixture = joinpath(@__DIR__, "..", "test", "fixtures", "spt3g_d1_full")
raw = readdlm(joinpath(fixture, "baseline_parameters.txt"), String; comments=true)
params_dict = Dict{String,Float64}(raw[i, 1] => parse(Float64, raw[i, 2]) for i in axes(raw, 1))
params = SPT3GD1Parameters(params_dict)
x0 = [params_dict[name] for name in SPT3G_D1_PARAMETER_NAMES]

cmb_mat = readdlm(joinpath(fixture, "cmb_spectra.txt"), Float64; comments=true)
cmb = SPT3GD1CMBTheory(Int.(cmb_mat[:, 1]), cmb_mat[:, 2], cmb_mat[:, 3], cmb_mat[:, 4])

like = SPT3GD1Likelihood()
model = SPT3GD1ForegroundModel()

_staged_predict(like, model, cmb, params) =
    bin_theory(like, instrument_stages(model, foreground_stages(model, cmb, params).dust_ee, params).calibration)

objective = function (x, fixed)
    l, m, c = fixed isa DifferentiationInterface.Constant ? fixed.data : fixed
    return loglikelihood(l, predict(l, m, c, SPT3GD1Parameters(x)))
end
fixed = DifferentiationInterface.Constant((like, model, cmb))

const N_ELL = length(model.ells)
const N_BANDPOWERS = length(like.data.data_vector)
const N_SPECTRA = length(like.data.windows)
println("Inputs  : $(N_SPECTRA) spectra, $(N_ELL) multipoles, $(N_BANDPOWERS) bandpowers, 43 nuisance parameters")

# ------------------------------------------------------ cold, fresh process ---
println("\n--- One-shot cold costs (fresh child process) ---")
cold = run_cold_harness()
for key in ("t_package_load_s", "t_artifact_construct_s", "t_first_predict_s",
            "t_prep_forwarddiff_s", "t_prep_mooncake_s")
    println(rpad(key, 26), round(cold[key], digits=3), " s")
end

# ---------------------------------------------------------- hot benchmarks ---
println("\n--- Hot BenchmarkTools measurements ---")
pred_baseline = predict(like, model, cmb, params)

b_predict_fast = @benchmark predict($like, $model, $cmb, $params)
b_predict_staged = @benchmark _staged_predict($like, $model, $cmb, $params)
b_loglike = @benchmark loglikelihood($like, $pred_baseline)
b_combined = @benchmark loglikelihood($like, predict($like, $model, $cmb, $params))

backend_mc = AutoMooncake(; config=nothing)
backend_fd = AutoForwardDiff()
prep_mc = DifferentiationInterface.prepare_gradient(objective, backend_mc, x0, fixed)
prep_fd = DifferentiationInterface.prepare_gradient(objective, backend_fd, x0, fixed)

# `evals=1`: the prepared caches are stateful and must not be reused within a sample.
b_mc_grad = @benchmark DifferentiationInterface.gradient($objective, $prep_mc, $backend_mc, $x0, $fixed) evals=1 seconds=15
b_fd_grad = @benchmark DifferentiationInterface.gradient($objective, $prep_fd, $backend_fd, $x0, $fixed) evals=1 seconds=15

const HOT = [
    ("Fast `predict` (production path)", b_predict_fast),
    ("Staged `predict` (readable path)", b_predict_staged),
    ("Data-only `loglikelihood`", b_loglike),
    ("Combined forward evaluation", b_combined),
    ("Hot prepared Mooncake gradient (43 params)", b_mc_grad),
    ("Hot prepared ForwardDiff gradient (43 params)", b_fd_grad),
]

for (name, b) in HOT
    println(rpad(name, 46), BenchmarkTools.prettytime(median(b).time))
end

# ------------------------------------------------------------- results.md ----
bytes_mib(x) = string(round(x / 1024^2, digits=2), " MiB")

open(RESULTS_PATH, "w") do io
    println(io, "<!-- Generated by benchmark/run_benchmarks.jl. Do not edit by hand. -->")
    println(io, "<!-- README.md and BENCHMARKS.md must quote this table verbatim. -->")
    println(io)
    println(io, "### Environment")
    println(io)
    println(io, "| Field | Value |")
    println(io, "|---|---|")
    println(io, "| Date | ", env.date, " |")
    println(io, "| Julia | ", env.julia, " |")
    println(io, "| OS | ", env.os, " |")
    println(io, "| CPU | ", env.cpu, " |")
    println(io, "| Logical cores | ", env.cores, " |")
    println(io, "| Julia threads | ", env.threads, " |")
    println(io, "| BLAS threads | ", env.blas_threads, " |")
    for (name, version) in sort(collect(env.versions))
        println(io, "| ", name, " | ", version, " |")
    end
    println(io, "| Inputs | ", N_SPECTRA, " spectra, ", N_ELL, " multipoles, ",
            N_BANDPOWERS, " bandpowers, 43 nuisance parameters |")
    println(io)
    println(io, "### Hot path (BenchmarkTools, compilation excluded)")
    println(io)
    println(io, "| Operation | Median | Minimum | Allocations | Memory |")
    println(io, "|---|---:|---:|---:|---:|")
    for (name, b) in HOT
        println(io, "| ", name,
                " | ", BenchmarkTools.prettytime(median(b).time),
                " | ", BenchmarkTools.prettytime(minimum(b).time),
                " | ", median(b).allocs,
                " | ", BenchmarkTools.prettymemory(median(b).memory), " |")
    end
    println(io)
    println(io, "### One-shot cold costs (fresh process, compilation included)")
    println(io)
    println(io, "These are single measurements from `benchmark/cold_prepare.jl`, timed with")
    println(io, "`time_ns`/`@timed` and peak RSS growth from `Sys.maxrss`. They are not")
    println(io, "BenchmarkTools statistics and have no median or minimum.")
    println(io)
    println(io, "| One-shot cost | Wall time | Allocated | Peak RSS growth |")
    println(io, "|---|---:|---:|---:|")
    println(io, "| Package load (`using SPTLikelihoods`) | ",
            round(cold["t_package_load_s"], digits=3), " s | — | ",
            bytes_mib(cold["rss_after_load_bytes"]), " |")
    println(io, "| Artifact-backed construction | ",
            round(cold["t_artifact_construct_s"], digits=3), " s | — | — |")
    println(io, "| First `predict` (compilation included) | ",
            round(cold["t_first_predict_s"], digits=3), " s | — | — |")
    println(io, "| ForwardDiff `prepare_gradient` | ",
            round(cold["t_prep_forwarddiff_s"], digits=3), " s | ",
            bytes_mib(cold["bytes_prep_forwarddiff"]), " | ",
            bytes_mib(cold["rss_prep_forwarddiff_bytes"]), " |")
    println(io, "| Mooncake `prepare_gradient` (tape build) | ",
            round(cold["t_prep_mooncake_s"], digits=3), " s | ",
            bytes_mib(cold["bytes_prep_mooncake"]), " | ",
            bytes_mib(cold["rss_prep_mooncake_bytes"]), " |")
    println(io, "| Whole cold process, load to prepared tape | ",
            round(cold["t_total_process_s"], digits=3), " s | — | ",
            bytes_mib(cold["rss_total_bytes"]), " |")
end

println("\nWrote authoritative results table to ", RESULTS_PATH)

# Keep README.md and BENCHMARKS.md quoting exactly this table.
include(joinpath(@__DIR__, "sync_docs.jl"))
for target in TARGETS
    println(sync_document(target) ? "updated  " : "unchanged", "  ", target)
end
