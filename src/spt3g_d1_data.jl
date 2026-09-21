const SPT3G_D1_SPECTRUM_ORDER = [
    "TT 90x90", "TE 90x90", "EE 90x90",
    "TT 90x150", "TE 90x150", "TE 150x90", "EE 90x150",
    "TT 90x220", "TE 90x220", "TE 220x90", "EE 90x220",
    "TT 150x150", "TE 150x150", "EE 150x150",
    "TT 150x220", "TE 150x220", "TE 220x150", "EE 150x220",
    "TT 220x220", "TE 220x220", "EE 220x220",
]

const SPT3G_D1_BINS_PER_SPECTRUM = [
    52, 72, 72, 52, 72, 72, 72, 52, 72, 72, 72,
    52, 72, 72, 52, 72, 72, 72, 52, 72, 72,
]

"""
    SPT3GD1Data

Fixed data and window functions for the full SPT-3G D1 T&E likelihood.
The unbinned theory convention is Dℓ on the explicitly stored multipole grid.
"""
struct SPT3GD1Data
    data_vector::Vector{Float64}
    cov_chol::Cholesky{Float64, Matrix{Float64}}
    windows::Vector{Matrix{Float64}}
    ells::Vector{Int}
    spectrum_order::Vector{String}
    bins_per_spectrum::Vector{Int}
end

function SPT3GD1Data(
    data_vector::AbstractVector{<:Real},
    covariance::AbstractMatrix{<:Real},
    windows::AbstractVector{<:AbstractMatrix{<:Real}},
    ells::AbstractVector{<:Integer};
    spectrum_order::AbstractVector{<:AbstractString}=SPT3G_D1_SPECTRUM_ORDER,
    bins_per_spectrum::AbstractVector{<:Integer}=SPT3G_D1_BINS_PER_SPECTRUM,
)
    n_spectra = length(spectrum_order)
    length(bins_per_spectrum) == n_spectra ||
        throw(DimensionMismatch("one bin count is required for each spectrum"))
    length(windows) == n_spectra ||
        throw(DimensionMismatch("one window matrix is required for each spectrum"))
    sum(bins_per_spectrum) == length(data_vector) ||
        throw(DimensionMismatch("data-vector length must equal the total number of bins"))
    size(covariance) == (length(data_vector), length(data_vector)) ||
        throw(DimensionMismatch("covariance dimensions must match the data vector"))
    issorted(ells) && all(diff(ells) .== 1) ||
        throw(ArgumentError("theory multipoles must be contiguous and sorted"))

    all(isfinite, data_vector) || throw(ArgumentError("data vector must be finite"))
    all(isfinite, covariance) || throw(ArgumentError("covariance must be finite"))

    n_ell = length(ells)
    typed_windows = Matrix{Float64}[]
    for (window, n_bins) in zip(windows, bins_per_spectrum)
        size(window) == (n_ell, n_bins) ||
            throw(DimensionMismatch("window dimensions must be (n_ell, n_bins)"))
        all(isfinite, window) || throw(ArgumentError("window matrices must be finite"))
        push!(typed_windows, Matrix{Float64}(window))
    end

    typed_covariance = Matrix{Float64}(covariance)
    isapprox(typed_covariance, typed_covariance'; rtol=0, atol=1e-12) ||
        throw(ArgumentError("covariance must be symmetric"))
    cov_chol = try
        cholesky(Symmetric(typed_covariance, :L))
    catch e
        e isa PosDefException && throw(ArgumentError("covariance must be positive definite"))
        rethrow(e)
    end

    return SPT3GD1Data(
        Vector{Float64}(data_vector),
        cov_chol,
        typed_windows,
        Vector{Int}(ells),
        String.(spectrum_order),
        Int.(bins_per_spectrum),
    )
end

"""Artifact name the converted D1 runtime tree must declare."""
const SPT3G_D1_ARTIFACT_NAME = "SPT3G_D1_TnE_v0"

"""Pinned `spt_candl_data` revision the D1 runtime artifact was converted from."""
const SPT3G_D1_SOURCE_REVISION = "bfe809a140087d19412aa6bc8c8d4ba18840b315"

"""Pinned `candl` revision of the engine this implementation reproduces."""
const SPT3G_D1_CANDL_REVISION = "650db0a6a0a2febed1e57350f0b991b537a4d2f7"

const SPT3G_D1_ELL_MIN = 2
const SPT3G_D1_ELL_MAX = 4095

_d1_meta_error(message) = throw(ArgumentError("SPT-3G D1 metadata.json: " * message))

function _json_capture(content::AbstractString, key::AbstractString, value_pattern::AbstractString)
    m = match(Regex("\"" * key * "\"\\s*:\\s*" * value_pattern), content)
    isnothing(m) && _d1_meta_error("missing or malformed required key \"$key\"")
    return m.captures[1]
end

_json_string(content, key) = String(_json_capture(content, key, "\"([^\"]*)\""))
_json_int(content, key) = parse(Int, _json_capture(content, key, "(-?\\d+)"))

function _json_string_array(content, key)
    body = _json_capture(content, key, "\\[([^\\]]*)\\]")
    return [String(strip(item, [' ', '\t', '\n', '\r', '"']))
            for item in split(body, ',') if !isempty(strip(item))]
end

function _json_int_array(content, key)
    body = _json_capture(content, key, "\\[([^\\]]*)\\]")
    return [parse(Int, strip(item)) for item in split(body, ',') if !isempty(strip(item))]
end

"""
    _validate_d1_metadata(meta_path) -> NamedTuple

Validate the converted D1 runtime tree's `metadata.json` against the release
constants this implementation was written for, and return the declared array
shapes so the caller can cross-check what it actually loaded.

Every check is mandatory: a key that is absent or malformed is an error rather
than a silently skipped validation. Validated here are the artifact identity,
the pinned `spt_candl_data` provenance revision and descriptor hash, the
multipole range, the 21-spectrum release ordering, the per-spectrum bin counts,
and the declared array dimensions.
"""
function _validate_d1_metadata(meta_path::AbstractString)
    isfile(meta_path) ||
        _d1_meta_error("no metadata.json at $meta_path; this is not a converted D1 artifact root")
    content = read(meta_path, String)

    # Artifact identity
    name = _json_string(content, "artifact")
    name == SPT3G_D1_ARTIFACT_NAME ||
        _d1_meta_error("declares artifact \"$name\", expected \"$SPT3G_D1_ARTIFACT_NAME\"")

    # Provenance: the pinned upstream revision and the likelihood descriptor hash
    revision = _json_string(content, "source_revision")
    revision == SPT3G_D1_SOURCE_REVISION ||
        _d1_meta_error("was converted from spt_candl_data $revision, expected $SPT3G_D1_SOURCE_REVISION")
    descriptor = _json_string(content, "descriptor_sha256")
    (length(descriptor) == 64 && all(isxdigit, descriptor)) ||
        _d1_meta_error("descriptor_sha256 \"$descriptor\" is not a SHA-256 digest")

    # Multipole support
    ell_min = _json_int(content, "ell_min")
    ell_max = _json_int(content, "ell_max")
    ell_min == SPT3G_D1_ELL_MIN ||
        _d1_meta_error("declares ell_min $ell_min, expected $SPT3G_D1_ELL_MIN")
    ell_max == SPT3G_D1_ELL_MAX ||
        _d1_meta_error("declares ell_max $ell_max, expected $SPT3G_D1_ELL_MAX")
    n_ell = ell_max - ell_min + 1

    # Release spectrum ordering and bin counts
    order = _json_string_array(content, "spectrum_order")
    order == SPT3G_D1_SPECTRUM_ORDER ||
        _d1_meta_error("declares a spectrum order that is not the 21-spectrum D1 release ordering")
    bins = _json_int_array(content, "bins_per_spectrum")
    bins == SPT3G_D1_BINS_PER_SPECTRUM ||
        _d1_meta_error("declares bins_per_spectrum $bins, expected $SPT3G_D1_BINS_PER_SPECTRUM")

    # Declared array dimensions
    n_bandpowers = sum(bins)
    declared_data = _json_int_array(content, "data_vector")
    declared_data == [n_bandpowers] ||
        _d1_meta_error("declares data_vector shape $declared_data, expected [$n_bandpowers]")
    declared_cov = _json_int_array(content, "covariance")
    declared_cov == [n_bandpowers, n_bandpowers] ||
        _d1_meta_error("declares covariance shape $declared_cov, expected [$n_bandpowers, $n_bandpowers]")
    # "windows" declares [n_spectra, n_ell, [bins...]], so only its leading two
    # integers are a flat array; read them directly.
    m_windows = match(r"\"windows\"\s*:\s*\[\s*(\d+)\s*,\s*(\d+)", content)
    isnothing(m_windows) && _d1_meta_error("missing or malformed required key \"windows\"")
    declared_windows = (parse(Int, m_windows.captures[1]), parse(Int, m_windows.captures[2]))
    declared_windows == (length(order), n_ell) ||
        _d1_meta_error("declares windows leading shape $declared_windows, expected $((length(order), n_ell))")

    return (
        artifact=name,
        source_revision=revision,
        descriptor_sha256=descriptor,
        ell_min=ell_min,
        ell_max=ell_max,
        n_ell=n_ell,
        spectrum_order=order,
        bins_per_spectrum=bins,
        n_bandpowers=n_bandpowers,
    )
end

"""
    _check_d1_against_metadata(meta, data_vector, covariance, windows, ells)

Cross-check arrays actually loaded from an artifact root against the shapes its
`metadata.json` declared. `SPT3GD1Data` separately enforces internal
consistency; this catches a tree whose arrays disagree with its own manifest.
"""
function _check_d1_against_metadata(meta, data_vector, covariance, windows, ells)
    isnothing(meta) && return nothing
    length(data_vector) == meta.n_bandpowers ||
        _d1_meta_error("data_vector has $(length(data_vector)) entries, metadata declares $(meta.n_bandpowers)")
    size(covariance) == (meta.n_bandpowers, meta.n_bandpowers) ||
        _d1_meta_error("covariance is $(size(covariance)), metadata declares $((meta.n_bandpowers, meta.n_bandpowers))")
    length(windows) == length(meta.spectrum_order) ||
        _d1_meta_error("$(length(windows)) window matrices loaded, metadata declares $(length(meta.spectrum_order))")
    length(ells) == meta.n_ell ||
        _d1_meta_error("$(length(ells)) multipoles loaded, metadata declares $(meta.n_ell)")
    (first(ells) == meta.ell_min && last(ells) == meta.ell_max) ||
        _d1_meta_error("multipole grid spans $(first(ells)):$(last(ells)), metadata declares $(meta.ell_min):$(meta.ell_max)")
    for (index, (window, n_bins)) in enumerate(zip(windows, meta.bins_per_spectrum))
        size(window) == (meta.n_ell, n_bins) ||
            _d1_meta_error("window $index is $(size(window)), metadata declares $((meta.n_ell, n_bins))")
    end
    return nothing
end

"""
    load_spt3g_d1_data(data_dir; require_metadata=false) -> SPT3GD1Data

Load the numeric artifact layout produced by
`validation/convert_spt3g_d1_artifact.py`. `data_dir` is the artifact root,
not the original candl release directory.

When `data_dir` contains a `metadata.json` it is validated by
[`_validate_d1_metadata`](@ref) and the loaded arrays are cross-checked against
the shapes it declares. Pass `require_metadata=true` to reject a tree that has
no manifest at all; the artifact-backed constructors do this, while the reduced
text fixtures used in the test suite do not carry one.
"""
function load_spt3g_d1_data(data_dir::AbstractString; require_metadata::Bool=false)
    meta_path = joinpath(data_dir, "metadata.json")
    meta = if isfile(meta_path)
        _validate_d1_metadata(meta_path)
    elseif require_metadata
        _d1_meta_error("no metadata.json at $meta_path; this is not a converted D1 artifact root")
    else
        nothing
    end
    if isfile(joinpath(data_dir, "data_vector.npy"))
        data_vector = vec(npzread(joinpath(data_dir, "data_vector.npy")))
        covariance = npzread(joinpath(data_dir, "covariance.npy"))
        ells = vec(npzread(joinpath(data_dir, "ells.npy")))
        windows = [
            npzread(joinpath(data_dir, "windows", "$(lpad(index, 2, '0')).npy"))
            for index in eachindex(SPT3G_D1_SPECTRUM_ORDER)
        ]
        _check_d1_against_metadata(meta, data_vector, covariance, windows, ells)
        return SPT3GD1Data(data_vector, covariance, windows, ells)
    elseif isfile(joinpath(data_dir, "data_vector.txt"))
        data_vector = vec(readdlm(joinpath(data_dir, "data_vector.txt"), Float64; comments=true, comment_char='#'))
        covariance = readdlm(joinpath(data_dir, "covariance.txt"), Float64; comments=true, comment_char='#')
        ells = Int.(vec(readdlm(joinpath(data_dir, "ells.txt"), Float64; comments=true, comment_char='#')))
        win_mat = readdlm(joinpath(data_dir, "windows.txt"), Float64; comments=true, comment_char='#')
        windows = [reshape(win_mat[:, index], :, 1) for index in 1:21]
        bins_per_spec = fill(1, 21)
        return SPT3GD1Data(data_vector, covariance, windows, ells; bins_per_spectrum=bins_per_spec)
    else
        throw(ArgumentError("No SPT-3G D1 data files found in $data_dir"))
    end
end

using Pkg.Artifacts: ensure_artifact_installed

function _get_d1_artifact_dir()
    artifacts_toml = joinpath(@__DIR__, "..", "Artifacts.toml")
    isfile(artifacts_toml) || return nothing
    h = Artifacts.artifact_hash("SPT3G_D1_TnE_v0", artifacts_toml)
    isnothing(h) && return nothing
    ensure_artifact_installed("SPT3G_D1_TnE_v0", artifacts_toml)
    return Artifacts.artifact_path(h)
end
