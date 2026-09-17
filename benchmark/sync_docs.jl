# Inject benchmark/results.md into the documents that quote it.
#
#     julia benchmark/sync_docs.jl
#
# `benchmark/results.md` is the single authoritative benchmark table. README.md
# and BENCHMARKS.md each carry a BEGIN/END marker pair whose contents are
# replaced from that file, so the two documents cannot drift apart or fall out of
# date relative to the last recorded run. `test/test_benchmark_docs.jl` asserts
# they are in sync, so hand-editing a number fails the suite.

const BEGIN_MARKER = "<!-- BEGIN BENCHMARK TABLE -->"
const END_MARKER = "<!-- END BENCHMARK TABLE -->"

const REPO_ROOT = normpath(joinpath(@__DIR__, ".."))
const RESULTS_PATH = joinpath(REPO_ROOT, "benchmark", "results.md")
const TARGETS = [joinpath(REPO_ROOT, "README.md"), joinpath(REPO_ROOT, "BENCHMARKS.md")]

"""
    benchmark_block() -> String

The authoritative table, wrapped in its markers.
"""
function benchmark_block()
    isfile(RESULTS_PATH) ||
        error("benchmark/results.md is missing; run `julia --project=benchmark benchmark/run_benchmarks.jl` first")
    return BEGIN_MARKER * "\n\n" * strip(read(RESULTS_PATH, String)) * "\n\n" * END_MARKER
end

"""
    sync_document(path) -> Bool

Replace the marked region of `path` with the authoritative table. Returns
whether the file changed.
"""
function sync_document(path::AbstractString)
    text = read(path, String)
    start_index = findfirst(BEGIN_MARKER, text)
    stop_index = findfirst(END_MARKER, text)
    (isnothing(start_index) || isnothing(stop_index)) &&
        error("$path has no $BEGIN_MARKER / $END_MARKER pair")
    updated = text[1:first(start_index)-1] * benchmark_block() * text[last(stop_index)+1:end]
    updated == text && return false
    write(path, updated)
    return true
end

if abspath(PROGRAM_FILE) == @__FILE__
    for target in TARGETS
        changed = sync_document(target)
        println(changed ? "updated  " : "unchanged", "  ", target)
    end
end
