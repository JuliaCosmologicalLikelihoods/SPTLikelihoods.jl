#!/usr/bin/env julia

using Pkg.Artifacts
using SHA

function usage()
    println("usage: julia --project=. validation/bind_spt3g_d1_artifact.jl <converted-data-dir> <archive-path> <zenodo-url>")
end

length(ARGS) == 3 || (usage(); exit(1))

source_dir = abspath(ARGS[1])
archive_path = abspath(ARGS[2])
zenodo_url = ARGS[3]
isdir(source_dir) || error("converted data directory not found: $source_dir")
isfile(archive_path) || error("artifact archive not found: $archive_path")

artifact_hash = create_artifact() do artifact_dir
    cp(source_dir, artifact_dir; force=true)
end
archive_sha256 = bytes2hex(open(sha256, archive_path))

bind_artifact!(
    joinpath(@__DIR__, "..", "Artifacts.toml"),
    "SPT3G_D1_TnE_v0",
    artifact_hash;
    download_info=[(zenodo_url, archive_sha256)],
    force=true,
)

println("artifact_git_tree_sha1=$(artifact_hash)")
println("archive_sha256=$archive_sha256")
