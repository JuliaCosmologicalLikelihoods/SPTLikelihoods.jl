const _D1_T2P_SIGMAS = (0.000274, 0.000192, 0.000169)

function _unpack_d1_spectra(long::AbstractVector, n_ell::Int)
    length(long) == 21 * n_ell ||
        throw(DimensionMismatch("D1 unbinned vector must contain 21 spectrum blocks"))
    T = eltype(long)
    TT = Array{T}(undef, 3, 3, n_ell)
    TE = Array{T}(undef, 3, 3, n_ell)
    EE = Array{T}(undef, 3, 3, n_ell)
    block(index) = @view long[(index - 1) * n_ell + 1:index * n_ell]

    for (index, (left, right)) in zip(_D1_TT_INDICES, _D1_TT_PAIRS)
        TT[left, right, :] = block(index)
        TT[right, left, :] = block(index)
    end
    for (index, (left, right)) in zip(_D1_TE_INDICES, _D1_TE_PAIRS)
        TE[left, right, :] = block(index)
    end
    for (index, (left, right)) in zip(_D1_EE_INDICES, _D1_EE_PAIRS)
        EE[left, right, :] = block(index)
        EE[right, left, :] = block(index)
    end
    return TT, TE, EE
end

function _pack_d1_spectra(TT::AbstractArray{T,3}, TE::AbstractArray{T,3}, EE::AbstractArray{T,3}) where T
    n_ell = size(TT, 3)
    size(TT) == (3, 3, n_ell) == size(TE) == size(EE) ||
        throw(DimensionMismatch("TT, TE, and EE tensors must each have shape (3, 3, n_ell)"))
    blocks = Vector{AbstractVector{T}}(undef, 21)
    for (index, (left, right)) in zip(_D1_TT_INDICES, _D1_TT_PAIRS)
        blocks[index] = @view TT[left, right, :]
    end
    for (index, (left, right)) in zip(_D1_TE_INDICES, _D1_TE_PAIRS)
        blocks[index] = @view TE[left, right, :]
    end
    for (index, (left, right)) in zip(_D1_EE_INDICES, _D1_EE_PAIRS)
        blocks[index] = @view EE[left, right, :]
    end
    return vcat(blocks...)
end

function _apply_d1_leakage(
    ells::AbstractVector,
    TT::AbstractArray{T,3},
    TE::AbstractArray{T,3},
    EE::AbstractArray{T,3},
    parameters::SPT3GD1Parameters,
) where T
    gamma = [parameters.t2p[index] .* (_D1_T2P_SIGMAS[index]^2 .* ells .^ 2) for index in 1:3]
    T_out = promote_type(T, typeof(parameters.t2p[1]))
    leaked_TE = similar(TE, T_out)
    leaked_EE = similar(EE, T_out)
    for right in 1:3, left in 1:3
        leaked_TE[left, right, :] = CMBForegrounds.te_leakage(
            @view(TT[left, right, :]), gamma[right],
        ) .+ @view(TE[left, right, :])
        leaked_EE[left, right, :] = CMBForegrounds.ee_leakage(
            @view(TT[left, right, :]), @view(TE[left, right, :]),
            @view(TE[right, left, :]), gamma[left], gamma[right],
        ) .+ @view(EE[left, right, :])
    end
    return TT, leaked_TE, leaked_EE
end

function _d1_temperature_beam(model::SPT3GD1ForegroundModel, parameters::SPT3GD1Parameters)
    n_ell = length(model.ells)
    T = promote_type(eltype(model.beam_modes), typeof(parameters.beam_modes[1]))
    response = ones(T, 3, n_ell)
    for mode in 1:9
        response .+= parameters.beam_modes[mode] .* @view(model.beam_modes[:, :, mode])
    end
    return response
end

function _d1_polarization_beam(model::SPT3GD1ForegroundModel, temperature_beam, parameters::SPT3GD1Parameters)
    n_ell = length(model.ells)
    T = promote_type(eltype(temperature_beam), typeof(parameters.beam_polarization[1]))
    response = Matrix{T}(undef, 3, n_ell)
    ell_800 = model.ell_800_index
    for frequency in 1:3
        beta = parameters.beam_polarization[frequency]
        main_beam = @view model.main_temperature_beams[frequency, :]
        normalization = main_beam[ell_800] + beta * (1 - main_beam[ell_800])
        response[frequency, :] = @. (main_beam + beta * (temperature_beam[frequency, :] - main_beam)) / normalization
    end
    return response
end

function _apply_d1_beams(
    TT::AbstractArray{T,3}, TE::AbstractArray{T,3}, EE::AbstractArray{T,3},
    temperature_beam::AbstractMatrix, polarization_beam::AbstractMatrix,
) where T
    T_out = promote_type(T, eltype(temperature_beam), eltype(polarization_beam))
    beamed_TT = similar(TT, T_out)
    beamed_TE = similar(TE, T_out)
    beamed_EE = similar(EE, T_out)
    for right in 1:3, left in 1:3
        beamed_TT[left, right, :] = @. TT[left, right, :] * temperature_beam[left, :] * temperature_beam[right, :]
        beamed_TE[left, right, :] = @. TE[left, right, :] * temperature_beam[left, :] * polarization_beam[right, :]
        beamed_EE[left, right, :] = @. EE[left, right, :] * polarization_beam[left, :] * polarization_beam[right, :]
    end
    return beamed_TT, beamed_TE, beamed_EE
end

function _apply_d1_calibration(
    TT::AbstractArray{T,3}, TE::AbstractArray{T,3}, EE::AbstractArray{T,3},
    parameters::SPT3GD1Parameters,
) where T
    temperature_calibration = (
        parameters.tcal_ext150 * parameters.tcal_rel90,
        parameters.tcal_ext150,
        parameters.tcal_ext150 * parameters.tcal_rel220,
    )
    polarization_calibration = (
        temperature_calibration[1] * parameters.ecal_ext150 * parameters.ecal_rel90,
        temperature_calibration[2] * parameters.ecal_ext150,
        temperature_calibration[3] * parameters.ecal_ext150 * parameters.ecal_rel220,
    )
    T_out = promote_type(T, typeof(parameters.tcal_ext150), typeof(parameters.ecal_ext150))
    calibrated_TT = similar(TT, T_out)
    calibrated_TE = similar(TE, T_out)
    calibrated_EE = similar(EE, T_out)
    for right in 1:3, left in 1:3
        calibrated_TT[left, right, :] = TT[left, right, :] ./ (temperature_calibration[left] * temperature_calibration[right])
        calibrated_TE[left, right, :] = TE[left, right, :] ./ (temperature_calibration[left] * polarization_calibration[right])
        calibrated_EE[left, right, :] = EE[left, right, :] ./ (polarization_calibration[left] * polarization_calibration[right])
    end
    return calibrated_TT, calibrated_TE, calibrated_EE
end

"""
    instrument_stages(model, preinstrument_Dls, parameters)

Apply the ordered D1 T-to-P leakage, polarized RC3 beam model, and map-leg
calibration to a release-ordered unbinned spectrum vector.
"""
function instrument_stages(
    model::SPT3GD1ForegroundModel,
    preinstrument_Dls::AbstractVector,
    parameters::SPT3GD1Parameters,
)
    TT, TE, EE = _unpack_d1_spectra(preinstrument_Dls, length(model.ells))
    leaked_TT, leaked_TE, leaked_EE = _apply_d1_leakage(model.ells, TT, TE, EE, parameters)
    leakage = _pack_d1_spectra(leaked_TT, leaked_TE, leaked_EE)
    temperature_beam = _d1_temperature_beam(model, parameters)
    polarization_beam = _d1_polarization_beam(model, temperature_beam, parameters)
    beamed_TT, beamed_TE, beamed_EE = _apply_d1_beams(leaked_TT, leaked_TE, leaked_EE, temperature_beam, polarization_beam)
    beam = _pack_d1_spectra(beamed_TT, beamed_TE, beamed_EE)
    calibrated_TT, calibrated_TE, calibrated_EE = _apply_d1_calibration(beamed_TT, beamed_TE, beamed_EE, parameters)
    calibration = _pack_d1_spectra(calibrated_TT, calibrated_TE, calibrated_EE)
    return (leakage=leakage, beam=beam, calibration=calibration)
end

function _instrument_parameter_values(parameters::SPT3GD1Parameters)
    return [
        parameters.t2p...,
        parameters.beam_modes...,
        parameters.beam_polarization...,
        parameters.tcal_ext150,
        parameters.tcal_rel90,
        parameters.tcal_rel220,
        parameters.ecal_ext150,
        parameters.ecal_rel90,
        parameters.ecal_rel220,
    ]
end

function _fast_instrument_core(
    model::SPT3GD1ForegroundModel,
    preinstrument_Dls::AbstractVector,
    parameters::AbstractVector,
)
    length(parameters) == 21 || throw(DimensionMismatch("D1 instrument parameter vector must have length 21"))
    TT, TE, EE = _unpack_d1_spectra(preinstrument_Dls, length(model.ells))
    gamma = [parameters[index] .* @view(model.t2p_kernels[index, :]) for index in 1:3]
    T = promote_type(eltype(preinstrument_Dls), eltype(parameters))
    temperature_beam = ones(T, 3, length(model.ells))
    for mode in 1:9
        temperature_beam .+= parameters[3 + mode] .* @view(model.beam_modes[:, :, mode])
    end
    polarization_beam = Matrix{T}(undef, 3, length(model.ells))
    for frequency in 1:3
        beta = parameters[12 + frequency]
        main_beam = @view model.main_temperature_beams[frequency, :]
        normalization = main_beam[model.ell_800_index] + beta * (1 - main_beam[model.ell_800_index])
        polarization_beam[frequency, :] = @. (main_beam + beta * (temperature_beam[frequency, :] - main_beam)) / normalization
    end

    temperature_calibration = (
        parameters[16] * parameters[17],
        parameters[16],
        parameters[16] * parameters[18],
    )
    polarization_calibration = (
        temperature_calibration[1] * parameters[19] * parameters[20],
        temperature_calibration[2] * parameters[19],
        temperature_calibration[3] * parameters[19] * parameters[21],
    )
    n_ell = length(model.ells)
    calibrated_TT = Array{T}(undef, 3, 3, n_ell)
    calibrated_TE = Array{T}(undef, 3, 3, n_ell)
    calibrated_EE = Array{T}(undef, 3, 3, n_ell)
    for right in 1:3, left in 1:3
        @. calibrated_TT[left, right, :] =
            TT[left, right, :] * temperature_beam[left, :] * temperature_beam[right, :] /
            (temperature_calibration[left] * temperature_calibration[right])
        @. calibrated_TE[left, right, :] =
            (TE[left, right, :] + gamma[right] * TT[left, right, :]) *
            temperature_beam[left, :] * polarization_beam[right, :] /
            (temperature_calibration[left] * polarization_calibration[right])
        @. calibrated_EE[left, right, :] =
            gamma[left] * TE[left, right, :] +
            gamma[right] * TE[right, left, :] +
            (gamma[left] * gamma[right]) * TT[left, right, :]
        @. calibrated_EE[left, right, :] =
            (EE[left, right, :] + calibrated_EE[left, right, :]) *
            polarization_beam[left, :] * polarization_beam[right, :] /
            (polarization_calibration[left] * polarization_calibration[right])
    end
    return _pack_d1_spectra(calibrated_TT, calibrated_TE, calibrated_EE)
end

function _fast_instrument_Dls(
    model::SPT3GD1ForegroundModel,
    preinstrument_Dls::AbstractVector,
    parameters::SPT3GD1Parameters,
)
    return _fast_instrument_core(model, preinstrument_Dls, _instrument_parameter_values(parameters))
end

function predict(
    like::SPT3GD1Likelihood,
    model::SPT3GD1ForegroundModel,
    cmb::SPT3GD1CMBTheory,
    parameters::SPT3GD1Parameters,
)
    _grids_match(cmb.ells, like.data.ells) || throw(ArgumentError("CMB ell grid must match likelihood ell grid"))
    _grids_match(model.ells, like.data.ells) || throw(ArgumentError("Foreground model ell grid must match likelihood ell grid"))
    preinstrument = _fast_preinstrument_Dls(model, cmb, parameters)
    return bin_theory(like, _fast_instrument_Dls(model, preinstrument, parameters))
end
