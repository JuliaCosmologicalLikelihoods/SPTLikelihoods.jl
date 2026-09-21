@testset "benchmark documentation is in sync" begin
    # `benchmark/results.md` is the single authoritative benchmark table.
    # README.md and BENCHMARKS.md must quote it verbatim, so a number edited by
    # hand into either document, or a stale table left behind after a rerun,
    # fails here rather than shipping two disagreeing tables.
    repo = normpath(joinpath(@__DIR__, ".."))
    results_path = joinpath(repo, "benchmark", "results.md")
    @test isfile(results_path)

    include(joinpath(repo, "benchmark", "sync_docs.jl"))
    block = benchmark_block()

    for name in ("README.md", "BENCHMARKS.md")
        text = read(joinpath(repo, name), String)
        @test occursin(BEGIN_MARKER, text)
        @test occursin(END_MARKER, text)
        @test occursin(block, text)
    end

    # The recorded table must describe the dependency this package actually
    # pins, not an older resolution.
    results = read(results_path, String)
    @test occursin("CMBForegrounds | 0.4", results)
end
