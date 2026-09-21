struct SPT3GD1CMBTheory{T}
    ells::Vector{Int}
    TT::Vector{T}
    TE::Vector{T}
    EE::Vector{T}
end

function SPT3GD1CMBTheory(
    ells::AbstractVector{<:Integer},
    TT::AbstractVector{T},
    TE::AbstractVector{T},
    EE::AbstractVector{T},
) where {T}
    n_ell = length(ells)
    length(TT) == n_ell == length(TE) == length(EE) ||
        throw(DimensionMismatch("TT, TE, EE, and ells must have the same length"))
    issorted(ells) && all(diff(ells) .== 1) ||
        throw(ArgumentError("theory multipoles must be contiguous and sorted"))
    all(isfinite, TT) && all(isfinite, TE) && all(isfinite, EE) ||
        throw(ArgumentError("theory spectra must be finite"))
    return SPT3GD1CMBTheory(Vector{Int}(ells), Vector{T}(TT), Vector{T}(TE), Vector{T}(EE))
end

struct SPT3GD1DustParameters{T}
    amplitude::T
    alpha::T
    beta::T
end

struct SPT3GD1Parameters{T}
    kappa::T
    poisson::NTuple{6,T}
    cib::NTuple{3,T}
    cib_alpha::T
    tsz::T
    ksz::T
    dust_tt::SPT3GD1DustParameters{T}
    dust_te::SPT3GD1DustParameters{T}
    dust_ee::SPT3GD1DustParameters{T}
    t2p::NTuple{3,T}
    beam_modes::NTuple{9,T}
    beam_polarization::NTuple{3,T}
    tcal_ext150::T
    tcal_rel90::T
    tcal_rel220::T
    ecal_ext150::T
    ecal_rel90::T
    ecal_rel220::T
end

function SPT3GD1Parameters(
    kappa,
    poisson::NTuple{6},
    cib::NTuple{3},
    cib_alpha,
    tsz,
    ksz,
    dust_tt::NTuple{3},
    dust_te::NTuple{3},
    dust_ee::NTuple{3},
    t2p::NTuple{3},
    beam_modes::NTuple{9},
    beam_polarization::NTuple{3},
    tcal_ext150,
    tcal_rel90,
    tcal_rel220,
    ecal_ext150,
    ecal_rel90,
    ecal_rel220,
)
    values = (
        kappa, poisson..., cib..., cib_alpha, tsz, ksz,
        dust_tt..., dust_te..., dust_ee..., t2p..., beam_modes...,
        beam_polarization..., tcal_ext150, tcal_rel90, tcal_rel220,
        ecal_ext150, ecal_rel90, ecal_rel220,
    )
    promoted = promote(values...)
    T = typeof(first(promoted))
    index = Ref(1)
    function next_parameter!()
        current = index[]
        index[] = current + 1
        return promoted[current]
    end
    kappa_t = next_parameter!()
    poisson_t = ntuple(_ -> next_parameter!(), 6)
    cib_t = ntuple(_ -> next_parameter!(), 3)
    cib_alpha_t, tsz_t, ksz_t = next_parameter!(), next_parameter!(), next_parameter!()
    dust_tt_t = SPT3GD1DustParameters(next_parameter!(), next_parameter!(), next_parameter!())
    dust_te_t = SPT3GD1DustParameters(next_parameter!(), next_parameter!(), next_parameter!())
    dust_ee_t = SPT3GD1DustParameters(next_parameter!(), next_parameter!(), next_parameter!())
    t2p_t = ntuple(_ -> next_parameter!(), 3)
    beam_modes_t = ntuple(_ -> next_parameter!(), 9)
    beam_polarization_t = ntuple(_ -> next_parameter!(), 3)
    return SPT3GD1Parameters{T}(
        kappa_t, poisson_t, cib_t, cib_alpha_t, tsz_t, ksz_t,
        dust_tt_t, dust_te_t, dust_ee_t, t2p_t, beam_modes_t,
        beam_polarization_t, next_parameter!(), next_parameter!(), next_parameter!(),
        next_parameter!(), next_parameter!(), next_parameter!(),
    )
end

_parameter_value(values::AbstractDict, name::AbstractString) = values[name]

function SPT3GD1Parameters(values::AbstractDict{<:AbstractString})
    poisson = ntuple(index -> _parameter_value(values, (
        "TT_Poisson_90x90", "TT_Poisson_90x150", "TT_Poisson_90x220",
        "TT_Poisson_150x150", "TT_Poisson_150x220", "TT_Poisson_220x220",
    )[index]), 6)
    cib = ntuple(index -> _parameter_value(values, (
        "TT_CIB_150x150", "TT_CIB_150x220", "TT_CIB_220x220",
    )[index]), 3)
    dust(prefix) = (
        _parameter_value(values, "$(prefix)_Amp"),
        _parameter_value(values, "$(prefix)_Alpha"),
        _parameter_value(values, "$(prefix)_Beta"),
    )
    return SPT3GD1Parameters(
        _parameter_value(values, "Kappa"), poisson, cib,
        _parameter_value(values, "TT_CIBClustering_Alpha"),
        _parameter_value(values, "TT_tSZ_Amp"), _parameter_value(values, "TT_kSZ_Amp"),
        dust("TT_GalCirrus"), dust("TE_PolGalDust"), dust("EE_PolGalDust"),
        ntuple(index -> _parameter_value(values, ("T2P2_90", "T2P2_150", "T2P2_220")[index]), 3),
        ntuple(index -> _parameter_value(values, "beta_$(index)"), 9),
        ntuple(index -> _parameter_value(values, ("beta_pol_90", "beta_pol_150", "beta_pol_220")[index]), 3),
        _parameter_value(values, "Tcal_ext150"),
        _parameter_value(values, "Tcal_rel90"),
        _parameter_value(values, "Tcal_rel220"),
        _parameter_value(values, "Ecal_ext150"),
        _parameter_value(values, "Ecal_rel90"),
        _parameter_value(values, "Ecal_rel220"),
    )
end

const SPT3G_D1_PARAMETER_NAMES = (
    "Kappa",
    "TT_Poisson_90x90", "TT_Poisson_90x150", "TT_Poisson_90x220",
    "TT_Poisson_150x150", "TT_Poisson_150x220", "TT_Poisson_220x220",
    "TT_CIB_150x150", "TT_CIB_150x220", "TT_CIB_220x220",
    "TT_CIBClustering_Alpha",
    "TT_tSZ_Amp", "TT_kSZ_Amp",
    "TT_GalCirrus_Amp", "TT_GalCirrus_Alpha", "TT_GalCirrus_Beta",
    "TE_PolGalDust_Amp", "TE_PolGalDust_Alpha", "TE_PolGalDust_Beta",
    "EE_PolGalDust_Amp", "EE_PolGalDust_Alpha", "EE_PolGalDust_Beta",
    "T2P2_90", "T2P2_150", "T2P2_220",
    "beta_1", "beta_2", "beta_3", "beta_4", "beta_5",
    "beta_6", "beta_7", "beta_8", "beta_9",
    "beta_pol_90", "beta_pol_150", "beta_pol_220",
    "Tcal_ext150", "Tcal_rel90", "Tcal_rel220",
    "Ecal_ext150", "Ecal_rel90", "Ecal_rel220",
)

const SPT3G_D1_FIDUCIAL_PARAMETERS = (
    Kappa = 0.0,
    TT_Poisson_90x90 = 8.0, TT_Poisson_90x150 = 10.0, TT_Poisson_90x220 = 17.0,
    TT_Poisson_150x150 = 15.0, TT_Poisson_150x220 = 25.0, TT_Poisson_220x220 = 60.0,
    TT_CIB_150x150 = 5.0, TT_CIB_150x220 = 15.0, TT_CIB_220x220 = 30.0,
    TT_CIBClustering_Alpha = 0.8,
    TT_tSZ_Amp = 4.0, TT_kSZ_Amp = 2.0,
    TT_GalCirrus_Amp = 1.5, TT_GalCirrus_Alpha = -2.5, TT_GalCirrus_Beta = 1.5,
    TE_PolGalDust_Amp = 0.1, TE_PolGalDust_Alpha = -2.4, TE_PolGalDust_Beta = 1.5,
    EE_PolGalDust_Amp = 0.05, EE_PolGalDust_Alpha = -2.4, EE_PolGalDust_Beta = 1.5,
    T2P2_90 = 0.0, T2P2_150 = 0.0, T2P2_220 = 0.0,
    beta_1 = 0.0, beta_2 = 0.0, beta_3 = 0.0, beta_4 = 0.0, beta_5 = 0.0,
    beta_6 = 0.0, beta_7 = 0.0, beta_8 = 0.0, beta_9 = 0.0,
    beta_pol_90 = 0.0, beta_pol_150 = 0.0, beta_pol_220 = 0.0,
    Tcal_ext150 = 1.0, Tcal_rel90 = 1.0, Tcal_rel220 = 1.0,
    Ecal_ext150 = 1.0, Ecal_rel90 = 1.0, Ecal_rel220 = 1.0,
)

"""
    SPT3GD1Parameters(values::NamedTuple)

Build the nuisance parameters from a name-keyed `NamedTuple` such as
[`SPT3G_D1_FIDUCIAL_PARAMETERS`](@ref). Keys follow
[`SPT3G_D1_PARAMETER_NAMES`](@ref).
"""
SPT3GD1Parameters(values::NamedTuple) =
    SPT3GD1Parameters(Dict{String,Float64}(string(k) => Float64(v) for (k, v) in pairs(values)))

"""
    SPT3GD1Parameters()

Return the fiducial 43-parameter nuisance set, i.e.
`SPT3GD1Parameters(SPT3G_D1_FIDUCIAL_PARAMETERS)`.
"""
SPT3GD1Parameters() = SPT3GD1Parameters(SPT3G_D1_FIDUCIAL_PARAMETERS)

function SPT3GD1Parameters(x::AbstractVector)
    length(x) == 43 || throw(DimensionMismatch("SPT-3G D1 parameter vector must have length 43, got $(length(x))"))
    return SPT3GD1Parameters(
        x[1],
        ntuple(i -> x[1 + i], 6),
        ntuple(i -> x[7 + i], 3),
        x[11], x[12], x[13],
        (x[14], x[15], x[16]),
        (x[17], x[18], x[19]),
        (x[20], x[21], x[22]),
        ntuple(i -> x[22 + i], 3),
        ntuple(i -> x[25 + i], 9),
        ntuple(i -> x[34 + i], 3),
        x[38], x[39], x[40], x[41], x[42], x[43],
    )
end

function _to_vector(p::SPT3GD1Parameters{T}) where T
    return T[
        p.kappa,
        p.poisson...,
        p.cib...,
        p.cib_alpha,
        p.tsz,
        p.ksz,
        p.dust_tt.amplitude, p.dust_tt.alpha, p.dust_tt.beta,
        p.dust_te.amplitude, p.dust_te.alpha, p.dust_te.beta,
        p.dust_ee.amplitude, p.dust_ee.alpha, p.dust_ee.beta,
        p.t2p...,
        p.beam_modes...,
        p.beam_polarization...,
        p.tcal_ext150, p.tcal_rel90, p.tcal_rel220,
        p.ecal_ext150, p.ecal_rel90, p.ecal_rel220,
    ]
end

Base.Vector(p::SPT3GD1Parameters) = _to_vector(p)
Base.collect(p::SPT3GD1Parameters) = _to_vector(p)

struct SPT3GD1ForegroundModel
    ells::Vector{Int}
    tsz_template::Vector{Float64}
    ksz_template::Vector{Float64}
    beam_modes::Array{Float64,3}
    main_temperature_beams::Matrix{Float64}
    ell_over_3000::Vector{Float64}
    ell_over_80::Vector{Float64}
    poisson_shape::Vector{Float64}
    t2p_kernels::Matrix{Float64}
    ell_800_index::Int
    tsz_scalings::NTuple{3,Float64}
end

function SPT3GD1ForegroundModel(
    ells::AbstractVector{<:Integer},
    tsz_template::AbstractVector{<:Real},
    ksz_template::AbstractVector{<:Real},
    ;
    beam_modes::AbstractArray{<:Real,3}=zeros(3, length(ells), 9),
    main_temperature_beams::AbstractMatrix{<:Real}=ones(3, length(ells)),
)
    length(ells) == length(tsz_template) == length(ksz_template) ||
        throw(DimensionMismatch("templates must share the theory ell grid"))
    pivot_index = findfirst(==(3000), ells)
    isnothing(pivot_index) && throw(ArgumentError("theory ell grid must include ell = 3000"))
    ell_800_index = findfirst(==(800), ells)
    isnothing(ell_800_index) && throw(ArgumentError("theory ell grid must include ell = 800"))
    size(beam_modes) == (3, length(ells), 9) ||
        throw(DimensionMismatch("beam modes must have shape (3, n_ell, 9)"))
    size(main_temperature_beams) == (3, length(ells)) ||
        throw(DimensionMismatch("main temperature beams must have shape (3, n_ell)"))
    all(isfinite, tsz_template) && all(isfinite, ksz_template) ||
        throw(ArgumentError("templates must be finite"))
    all(isfinite, beam_modes) && all(isfinite, main_temperature_beams) ||
        throw(ArgumentError("beam arrays must be finite"))
    tsz_norm = Float64.(tsz_template) ./ tsz_template[pivot_index]
    ksz_norm = Float64.(ksz_template) ./ ksz_template[pivot_index]
    ell_over_3000 = Float64.(ells) ./ 3000
    ell_over_80 = Float64.(ells) ./ 80
    poisson_shape = angular_power(PowerLawShape(3000.0), Float64.(ells), 2.0)
    t2p_kernels = stack(
        (_D1_T2P_SIGMAS[index]^2 .* Float64.(ells) .^ 2 for index in 1:3);
        dims=1,
    )
    tsz_scalings = ntuple(index -> _candl_tsz_scaling(_D1_TSZ_FREQUENCIES[index], 143.0), 3)
    return SPT3GD1ForegroundModel(
        Vector{Int}(ells), tsz_norm, ksz_norm, Array{Float64,3}(beam_modes),
        Matrix{Float64}(main_temperature_beams), ell_over_3000, ell_over_80,
        poisson_shape, t2p_kernels, ell_800_index, tsz_scalings,
    )
end

function SPT3GD1ForegroundModel(data_dir::AbstractString)
    if isfile(joinpath(data_dir, "ells.npy"))
        ells = vec(npzread(joinpath(data_dir, "ells.npy")))
        tsz_template = vec(npzread(joinpath(data_dir, "templates", "tsz.npy")))
        ksz_template = vec(npzread(joinpath(data_dir, "templates", "ksz.npy")))
        beam_modes = npzread(joinpath(data_dir, "beams", "eigenmodes.npy"))
        main_beams = npzread(joinpath(data_dir, "beams", "main_temperature.npy"))
        return SPT3GD1ForegroundModel(
            ells, tsz_template, ksz_template;
            beam_modes=beam_modes, main_temperature_beams=main_beams,
        )
    elseif isfile(joinpath(data_dir, "cmb_spectra.txt"))
        cmb_mat = readdlm(joinpath(data_dir, "cmb_spectra.txt"), Float64; comments=true, comment_char='#')
        ells = Int.(cmb_mat[:, 1])
        fg_mat = readdlm(joinpath(data_dir, "foreground_templates.txt"), Float64; comments=true, comment_char='#')
        tsz_template = fg_mat[:, 1]
        ksz_template = fg_mat[:, 2]
        main_beams = Matrix(readdlm(joinpath(data_dir, "beam_main.txt"), Float64; comments=true, comment_char='#')')
        modes_mat = readdlm(joinpath(data_dir, "beam_modes.txt"), Float64; comments=true, comment_char='#')
        n_ell = length(ells)
        beam_modes = permutedims(reshape(modes_mat, n_ell, 9, 3), (3, 1, 2))
        return SPT3GD1ForegroundModel(
            ells, tsz_template, ksz_template;
            beam_modes=beam_modes, main_temperature_beams=main_beams,
        )
    else
        throw(ArgumentError("No SPT-3G D1 model files found in $data_dir"))
    end
end

function SPT3GD1ForegroundModel()
    artifact_dir = _get_d1_artifact_dir()
    if isnothing(artifact_dir) || !isdir(artifact_dir)
        throw(ArgumentError("SPT-3G D1 artifact \"SPT3G_D1_TnE_v0\" is not bound or installed in Artifacts.toml. Please provide the artifact path via SPT3GD1ForegroundModel(data_dir)."))
    end
    _validate_d1_metadata(joinpath(artifact_dir, "metadata.json"))
    return SPT3GD1ForegroundModel(artifact_dir)
end

const _D1_TT_INDICES = (1, 4, 8, 12, 15, 19)
const _D1_TE_INDICES = (2, 5, 6, 9, 10, 13, 16, 17, 20)
const _D1_EE_INDICES = (3, 7, 11, 14, 18, 21)
const _D1_TT_PAIRS = ((1, 1), (1, 2), (1, 3), (2, 2), (2, 3), (3, 3))
const _D1_TE_PAIRS = ((1, 1), (1, 2), (2, 1), (1, 3), (3, 1), (2, 2), (2, 3), (3, 2), (3, 3))
const _D1_EE_PAIRS = ((1, 1), (1, 2), (1, 3), (2, 2), (2, 3), (3, 3))
const _D1_TSZ_FREQUENCIES = (95.6933, 148.849, 220.15)
const _D1_CIRRUS_FREQUENCIES = (95.9515, 149.993, 222.745)
const _D1_POLARISED_DUST_FREQUENCIES = (95.9631, 150.012, 222.773)
const _D1_T_CMB = 2.72548
const _D1_GHZ_KELVIN = 6.62606957e-34 / 1.3806488e-23 * 1e9

function _candl_dcl_dell(ells::AbstractVector, Dls::AbstractVector)
    Cls = @. Dls * (2π) / (ells * (ells + 1))
    derivative = similar(Cls)
    derivative[1] = (Cls[2] - Cls[1]) / (ells[2] - ells[1])
    @inbounds for index in 2:length(ells)-1
        derivative[index] = (Cls[index + 1] - Cls[index - 1]) / (ells[index + 1] - ells[index - 1])
    end
    derivative[end] = (Cls[end] - Cls[end - 1]) / (ells[end] - ells[end - 1])
    return derivative
end

function _candl_apply_ssl(ells::AbstractVector, kappa, Dls::AbstractVector)
    derivative = _candl_dcl_dell(ells, Dls)
    return @. Dls - kappa * (ells * derivative * (ells * (ells + 1)) / (2π) + 2 * Dls)
end

function _candl_apply_aberration(ells::AbstractVector, coefficient, Dls::AbstractVector)
    derivative = _candl_dcl_dell(ells, Dls)
    return @. Dls - coefficient * ells * derivative * (ells * (ells + 1)) / (2π)
end

function _candl_tsz_scaling(nu, nu0)
    x0 = _D1_GHZ_KELVIN * nu0 / _D1_T_CMB
    x = _D1_GHZ_KELVIN * nu / _D1_T_CMB
    g0 = x0 * (exp(x0) + 1) / (exp(x0) - 1) - 4
    g = x * (exp(x) + 1) / (exp(x) - 1) - 4
    return g / g0
end

function _candl_blackbody(nu, nu0, temperature)
    x0 = _D1_GHZ_KELVIN * nu0 / temperature
    x = _D1_GHZ_KELVIN * nu / temperature
    return (nu / nu0)^3 * (exp(x0) - 1) / (exp(x) - 1)
end

function _candl_blackbody_derivative(nu, nu0, temperature)
    x0 = _D1_GHZ_KELVIN * nu0 / temperature
    x = _D1_GHZ_KELVIN * nu / temperature
    d0 = x0^4 * exp(x0) / (exp(x0) - 1)^2
    d = x^4 * exp(x) / (exp(x) - 1)^2
    return d / d0
end

function _candl_dust_power(ells::AbstractVector, dust::SPT3GD1DustParameters, nu1, nu2)
    sed(nu) = (nu / 150.0)^dust.beta *
              _candl_blackbody(nu, 150.0, 19.6) /
              _candl_blackbody_derivative(nu, 150.0, _D1_T_CMB)
    return @. dust.amplitude * sed(nu1) * sed(nu2) * (ells / 80)^((dust.alpha) + 2)
end

function _long_from_type_spectra(TT::AbstractVector, TE::AbstractVector, EE::AbstractVector)
    return vcat(TT, TE, EE, TT, TE, TE, EE, TT, TE, TE, EE, TT, TE, EE, TT, TE, TE, EE, TT, TE, EE)
end

function _empty_blocks(model::SPT3GD1ForegroundModel, parameter_type::Type)
    return [zeros(parameter_type, length(model.ells)) for _ in eachindex(SPT3G_D1_SPECTRUM_ORDER)]
end

function _blocks_long(blocks::AbstractVector)
    return vcat(blocks...)
end

function _tt_component(model::SPT3GD1ForegroundModel, parameter_type::Type, values::NTuple{6})
    blocks = _empty_blocks(model, parameter_type)
    for (index, amplitude) in zip(_D1_TT_INDICES, values)
        blocks[index] = amplitude
    end
    return blocks
end

function _poisson_long(model::SPT3GD1ForegroundModel, parameters::SPT3GD1Parameters{T}) where T
    blocks = _empty_blocks(model, T)
    for (index, amplitude) in zip(_D1_TT_INDICES, parameters.poisson)
        blocks[index] = amplitude .* model.poisson_shape
    end
    return _blocks_long(blocks)
end

function _cib_long(model::SPT3GD1ForegroundModel, parameters::SPT3GD1Parameters{T}) where T
    shape = angular_power(PowerLawShape(3000.0), model.ells, parameters.cib_alpha)
    blocks = _empty_blocks(model, T)
    for (index, amplitude) in zip((12, 15, 19), parameters.cib)
        blocks[index] = amplitude .* shape
    end
    return _blocks_long(blocks)
end

function _tsz_auto(model::SPT3GD1ForegroundModel, parameters::SPT3GD1Parameters, frequency::Real)
    scaling = _candl_tsz_scaling(frequency, 143.0)
    return @. parameters.tsz * scaling * scaling * model.tsz_template
end

function _tsz_long(model::SPT3GD1ForegroundModel, parameters::SPT3GD1Parameters{T}) where T
    blocks = _empty_blocks(model, T)
    for (index, (left, right)) in zip(_D1_TT_INDICES, _D1_TT_PAIRS)
        scaling = _candl_tsz_scaling(_D1_TSZ_FREQUENCIES[left], 143.0) *
                  _candl_tsz_scaling(_D1_TSZ_FREQUENCIES[right], 143.0)
        blocks[index] = @. parameters.tsz * scaling * model.tsz_template
    end
    return _blocks_long(blocks)
end

function _tsz_cib_long(model::SPT3GD1ForegroundModel, parameters::SPT3GD1Parameters{T}) where T
    cib_amplitudes = (zero(parameters.cib[1]), parameters.cib[1], parameters.cib[3])
    cib_shape = (model.ells ./ 3000) .^ parameters.cib_alpha
    tsz_auto = [_tsz_auto(model, parameters, frequency) for frequency in _D1_TSZ_FREQUENCIES]
    blocks = _empty_blocks(model, T)
    for (index, (left, right)) in zip(_D1_TT_INDICES, _D1_TT_PAIRS)
        term = zeros(T, length(model.ells))
        if left != 1
            cib_left = cib_amplitudes[left] .* cib_shape
            @. term += sqrt(cib_left * tsz_auto[right])
        end
        if right != 1
            cib_right = cib_amplitudes[right] .* cib_shape
            @. term += sqrt(cib_right * tsz_auto[left])
        end
        blocks[index] = @. -0.26 * term
    end
    return _blocks_long(blocks)
end

function _ksz_long(model::SPT3GD1ForegroundModel, parameters::SPT3GD1Parameters{T}) where T
    blocks = _empty_blocks(model, T)
    scaled = CMBForegrounds.ksz_template_scaled(model.ksz_template, parameters.ksz)
    for index in _D1_TT_INDICES
        blocks[index] = scaled
    end
    return _blocks_long(blocks)
end

function _dust_long(
    model::SPT3GD1ForegroundModel,
    parameters::SPT3GD1Parameters{T},
    dust::SPT3GD1DustParameters,
    frequencies::NTuple{3,Float64},
    indices,
    pairs,
) where T
    blocks = _empty_blocks(model, T)
    for (index, (left, right)) in zip(indices, pairs)
        blocks[index] = _candl_dust_power(model.ells, dust, frequencies[left], frequencies[right])
    end
    return _blocks_long(blocks)
end

"""
    foreground_stages(model, cmb, parameters)

Return the release-ordered unbinned spectra after every D1 transformation through
the three Galactic-dust stages. Leakage, beam modes, and calibration are applied
in the subsequent instrument layer.
"""
function foreground_stages(
    model::SPT3GD1ForegroundModel,
    cmb::SPT3GD1CMBTheory,
    parameters::SPT3GD1Parameters,
)
    _grids_match(cmb.ells, model.ells) || throw(ArgumentError("CMB and foreground-model ell grids differ"))
    input = _long_from_type_spectra(cmb.TT, cmb.TE, cmb.EE)
    ssl_tt = _candl_apply_ssl(model.ells, parameters.kappa, cmb.TT)
    ssl_te = _candl_apply_ssl(model.ells, parameters.kappa, cmb.TE)
    ssl_ee = _candl_apply_ssl(model.ells, parameters.kappa, cmb.EE)
    ssl = _long_from_type_spectra(ssl_tt, ssl_te, ssl_ee)
    aberration_tt = _candl_apply_aberration(model.ells, -0.0004826, ssl_tt)
    aberration_te = _candl_apply_aberration(model.ells, -0.0004826, ssl_te)
    aberration_ee = _candl_apply_aberration(model.ells, -0.0004826, ssl_ee)
    aberration = _long_from_type_spectra(aberration_tt, aberration_te, aberration_ee)
    poisson = aberration .+ _poisson_long(model, parameters)
    cib = poisson .+ _cib_long(model, parameters)
    tsz = cib .+ _tsz_long(model, parameters)
    tsz_cib = tsz .+ _tsz_cib_long(model, parameters)
    ksz = tsz_cib .+ _ksz_long(model, parameters)
    dust_tt = ksz .+ _dust_long(model, parameters, parameters.dust_tt, _D1_CIRRUS_FREQUENCIES, _D1_TT_INDICES, _D1_TT_PAIRS)
    dust_te = dust_tt .+ _dust_long(model, parameters, parameters.dust_te, _D1_POLARISED_DUST_FREQUENCIES, _D1_TE_INDICES, _D1_TE_PAIRS)
    dust_ee = dust_te .+ _dust_long(model, parameters, parameters.dust_ee, _D1_POLARISED_DUST_FREQUENCIES, _D1_EE_INDICES, _D1_EE_PAIRS)
    return (
        input=input,
        ssl=ssl,
        aberration=aberration,
        poisson=poisson,
        cib=cib,
        tsz=tsz,
        tsz_cib=tsz_cib,
        ksz=ksz,
        dust_tt=dust_tt,
        dust_te=dust_te,
        dust_ee=dust_ee,
    )
end

function _fast_preinstrument_Dls(
    model::SPT3GD1ForegroundModel,
    cmb::SPT3GD1CMBTheory,
    parameters::SPT3GD1Parameters,
)
    ssl_tt = _candl_apply_ssl(model.ells, parameters.kappa, cmb.TT)
    ssl_te = _candl_apply_ssl(model.ells, parameters.kappa, cmb.TE)
    ssl_ee = _candl_apply_ssl(model.ells, parameters.kappa, cmb.EE)
    aberration_tt = _candl_apply_aberration(model.ells, -0.0004826, ssl_tt)
    aberration_te = _candl_apply_aberration(model.ells, -0.0004826, ssl_te)
    aberration_ee = _candl_apply_aberration(model.ells, -0.0004826, ssl_ee)
    T = promote_type(eltype(cmb.TT), eltype(cmb.TE), eltype(cmb.EE), typeof(parameters.kappa), typeof(parameters.poisson[1]))
    output = Vector{T}(_long_from_type_spectra(aberration_tt, aberration_te, aberration_ee))
    ell = model.ells
    n_ell = length(ell)
    blocks = ntuple(index -> @view(output[(index - 1) * n_ell + 1:index * n_ell]), 21)
    poisson_shape = model.poisson_shape
    for (index, amplitude) in zip(_D1_TT_INDICES, parameters.poisson)
        @. blocks[index] += amplitude * poisson_shape
    end

    cib_shape = angular_power(PowerLawShape(3000.0), model.ells, parameters.cib_alpha)
    for (index, amplitude) in zip((12, 15, 19), parameters.cib)
        @. blocks[index] += amplitude * cib_shape
    end

    for (index, (left, right)) in zip(_D1_TT_INDICES, _D1_TT_PAIRS)
        scaling = model.tsz_scalings[left] * model.tsz_scalings[right]
        @. blocks[index] += parameters.tsz * scaling * model.tsz_template
    end

    cib_amplitudes = (zero(parameters.cib[1]), parameters.cib[1], parameters.cib[3])
    tsz_scalings = model.tsz_scalings
    for (index, (left, right)) in zip(_D1_TT_INDICES, _D1_TT_PAIRS)
        @inbounds for ell_index in eachindex(ell)
            correlation = zero(parameters.tsz)
            if left != 1
                cib_left = cib_amplitudes[left] * cib_shape[ell_index]
                tsz_right = parameters.tsz * tsz_scalings[right]^2 * model.tsz_template[ell_index]
                correlation += sqrt(cib_left * tsz_right)
            end
            if right != 1
                cib_right = cib_amplitudes[right] * cib_shape[ell_index]
                tsz_left = parameters.tsz * tsz_scalings[left]^2 * model.tsz_template[ell_index]
                correlation += sqrt(cib_right * tsz_left)
            end
            blocks[index][ell_index] += -0.26 * correlation
        end
    end

    for index in _D1_TT_INDICES
        @. blocks[index] += parameters.ksz * model.ksz_template
    end

    for (dust, frequencies, indices, pairs) in (
        (parameters.dust_tt, _D1_CIRRUS_FREQUENCIES, _D1_TT_INDICES, _D1_TT_PAIRS),
        (parameters.dust_te, _D1_POLARISED_DUST_FREQUENCIES, _D1_TE_INDICES, _D1_TE_PAIRS),
        (parameters.dust_ee, _D1_POLARISED_DUST_FREQUENCIES, _D1_EE_INDICES, _D1_EE_PAIRS),
    )
        dust_shape = model.ell_over_80 .^ (dust.alpha + 2)
        for (index, (left, right)) in zip(indices, pairs)
            sed_left = (frequencies[left] / 150.0)^dust.beta *
                       _candl_blackbody(frequencies[left], 150.0, 19.6) /
                       _candl_blackbody_derivative(frequencies[left], 150.0, _D1_T_CMB)
            sed_right = (frequencies[right] / 150.0)^dust.beta *
                        _candl_blackbody(frequencies[right], 150.0, 19.6) /
                        _candl_blackbody_derivative(frequencies[right], 150.0, _D1_T_CMB)
            @. blocks[index] += dust.amplitude * sed_left * sed_right * dust_shape
        end
    end
    return output
end
