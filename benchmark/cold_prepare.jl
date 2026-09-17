# Fresh-process cold-cost harness for the SPT-3G D1 benchmarks.
#
# Started as a child process by `benchmark/run_benchmarks.jl`. Everything here is
# measured exactly once, in a process that has never loaded the package, because
# these costs are one-shot by nature and cannot be measured by BenchmarkTools:
# repeating them in a warm process would measure something else entirely.
#
# Wall time comes from `time_ns()` and resident-memory growth from `Sys.maxrss()`,
# which reports the process peak RSS. Reported RSS figures are therefore the peak
# growth across the call, not an instantaneous reading. These are explicitly NOT
# BenchmarkTools statistics and are reported separately from the hot table.
#
# Emits `key=value` lines on stdout for the parent to parse.

const T_START = time_ns()
const RSS_START = Sys.maxrss()

using DelimitedFiles
t0 = time_ns()
using SPTLikelihoods
t_load = (time_ns() - t0) / 1e9
rss_after_load = Sys.maxrss()

using DifferentiationInterface
using ADTypes: AutoForwardDiff, AutoMooncake
import Mooncake
import ForwardDiff

# Artifact-backed construction, cold
t0 = time_ns()
like = SPT3GD1Likelihood()
model = SPT3GD1ForegroundModel()
t_construct = (time_ns() - t0) / 1e9

fixture = joinpath(@__DIR__, "..", "test", "fixtures", "spt3g_d1_full")
raw = readdlm(joinpath(fixture, "baseline_parameters.txt"), String; comments=true)
params_dict = Dict{String,Float64}(raw[i, 1] => parse(Float64, raw[i, 2]) for i in axes(raw, 1))
params = SPT3GD1Parameters(params_dict)
x0 = [params_dict[name] for name in SPT3G_D1_PARAMETER_NAMES]
cmb_mat = readdlm(joinpath(fixture, "cmb_spectra.txt"), Float64; comments=true)
cmb = SPT3GD1CMBTheory(Int.(cmb_mat[:, 1]), cmb_mat[:, 2], cmb_mat[:, 3], cmb_mat[:, 4])

# First `predict`: includes method compilation
t0 = time_ns()
pred = predict(like, model, cmb, params)
t_first_predict = (time_ns() - t0) / 1e9

objective = function (x, fixed)
    l, m, c = fixed isa DifferentiationInterface.Constant ? fixed.data : fixed
    return loglikelihood(l, predict(l, m, c, SPT3GD1Parameters(x)))
end
fixed = DifferentiationInterface.Constant((like, model, cmb))

# ForwardDiff preparation
rss_before_fd = Sys.maxrss()
stats_fd = @timed DifferentiationInterface.prepare_gradient(objective, AutoForwardDiff(), x0, fixed)
rss_after_fd = Sys.maxrss()

# Mooncake tape construction: the dominant one-time cost
rss_before_mc = Sys.maxrss()
stats_mc = @timed DifferentiationInterface.prepare_gradient(objective, AutoMooncake(; config=nothing), x0, fixed)
rss_after_mc = Sys.maxrss()

println("t_package_load_s=", t_load)
println("t_artifact_construct_s=", t_construct)
println("t_first_predict_s=", t_first_predict)
println("rss_after_load_bytes=", rss_after_load - RSS_START)
println("t_prep_forwarddiff_s=", stats_fd.time)
println("bytes_prep_forwarddiff=", stats_fd.bytes)
println("rss_prep_forwarddiff_bytes=", rss_after_fd - rss_before_fd)
println("t_prep_mooncake_s=", stats_mc.time)
println("bytes_prep_mooncake=", stats_mc.bytes)
println("rss_prep_mooncake_bytes=", rss_after_mc - rss_before_mc)
println("rss_total_bytes=", Sys.maxrss() - RSS_START)
println("t_total_process_s=", (time_ns() - T_START) / 1e9)
println("n_bandpowers=", length(pred))
