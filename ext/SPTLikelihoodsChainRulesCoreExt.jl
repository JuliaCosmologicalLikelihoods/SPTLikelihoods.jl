module SPTLikelihoodsChainRulesCoreExt

using SPTLikelihoods: _fixed_covariance_solve, _fast_instrument_core, _grids_match,
    _unpack_d1_spectra, _pack_d1_spectra, SPT3GD1ForegroundModel,
    _D1_TT_INDICES, _D1_TT_PAIRS, _D1_TE_INDICES, _D1_TE_PAIRS,
    _D1_EE_INDICES, _D1_EE_PAIRS
using LinearAlgebra: Cholesky
import ChainRulesCore
using ChainRulesCore: NoTangent, ProjectTo, unthunk

function ChainRulesCore.rrule(
    ::typeof(_fixed_covariance_solve),
    cov_chol::Cholesky,
    residual::AbstractVector,
)
    solved = _fixed_covariance_solve(cov_chol, residual)
    project_residual = ProjectTo(residual)

    function fixed_covariance_solve_pullback(solved_bar_thunked)
        solved_bar = unthunk(solved_bar_thunked)
        residual_bar = project_residual(cov_chol \ solved_bar)
        return NoTangent(), NoTangent(), residual_bar
    end

    return solved, fixed_covariance_solve_pullback
end

function ChainRulesCore.rrule(
    ::typeof(_fast_instrument_core),
    model::SPT3GD1ForegroundModel,
    preinstrument::AbstractVector,
    parameters::AbstractVector,
)
    output = _fast_instrument_core(model, preinstrument, parameters)
    TT, TE, EE = _unpack_d1_spectra(preinstrument, length(model.ells))
    n_ell = length(model.ells)
    gamma = [parameters[index] .* @view(model.t2p_kernels[index, :]) for index in 1:3]

    temperature_beam = ones(eltype(parameters), 3, n_ell)
    for mode in 1:9
        temperature_beam .+= parameters[3 + mode] .* @view(model.beam_modes[:, :, mode])
    end
    polarization_beam = similar(temperature_beam)
    for frequency in 1:3
        beta = parameters[12 + frequency]
        main_beam = @view model.main_temperature_beams[frequency, :]
        normalization = main_beam[model.ell_800_index] + beta * (1 - main_beam[model.ell_800_index])
        polarization_beam[frequency, :] = @. (main_beam + beta * (temperature_beam[frequency, :] - main_beam)) / normalization
    end
    temperature_calibration = (
        parameters[16] * parameters[17], parameters[16], parameters[16] * parameters[18],
    )
    polarization_calibration = (
        temperature_calibration[1] * parameters[19] * parameters[20],
        temperature_calibration[2] * parameters[19],
        temperature_calibration[3] * parameters[19] * parameters[21],
    )
    project_preinstrument = ProjectTo(preinstrument)
    project_parameters = ProjectTo(parameters)

    function fast_instrument_pullback(output_bar_thunked)
        output_bar = unthunk(output_bar_thunked)
        TT_bar, TE_bar, EE_bar = _unpack_d1_spectra(output_bar, n_ell)
        dTT = zeros(eltype(output_bar), 3, 3, n_ell)
        dTE = zeros(eltype(output_bar), 3, 3, n_ell)
        dEE = zeros(eltype(output_bar), 3, 3, n_ell)
        dgamma = [zeros(eltype(output_bar), n_ell) for _ in 1:3]
        dtemperature_beam = zeros(eltype(output_bar), 3, n_ell)
        dpolarization_beam = zeros(eltype(output_bar), 3, n_ell)
        dtemperature_calibration = zeros(eltype(output_bar), 3)
        dpolarization_calibration = zeros(eltype(output_bar), 3)

        for (index, (left, right)) in zip(_D1_TT_INDICES, _D1_TT_PAIRS)
            bar = @view TT_bar[left, right, :]
            denom = temperature_calibration[left] * temperature_calibration[right]
            factor = @. temperature_beam[left, :] * temperature_beam[right, :] / denom
            spectrum = @view TT[left, right, :]
            @. dTT[left, right, :] += bar * factor
            @. dtemperature_beam[left, :] += bar * spectrum * temperature_beam[right, :] / denom
            @. dtemperature_beam[right, :] += bar * spectrum * temperature_beam[left, :] / denom
            output_block = @. spectrum * factor
            dtemperature_calibration[left] -= sum(bar .* output_block) / temperature_calibration[left]
            dtemperature_calibration[right] -= sum(bar .* output_block) / temperature_calibration[right]
        end

        for (index, (left, right)) in zip(_D1_TE_INDICES, _D1_TE_PAIRS)
            bar = @view TE_bar[left, right, :]
            tt_left, tt_right = minmax(left, right)
            denom = temperature_calibration[left] * polarization_calibration[right]
            factor = @. temperature_beam[left, :] * polarization_beam[right, :] / denom
            TT_block = @view TT[tt_left, tt_right, :]
            TE_block = @view TE[left, right, :]
            leaked = @. TE_block + gamma[right] * TT_block
            @. dTE[left, right, :] += bar * factor
            @. dTT[tt_left, tt_right, :] += bar * factor * gamma[right]
            @. dgamma[right] += bar * factor * TT_block
            @. dtemperature_beam[left, :] += bar * leaked * polarization_beam[right, :] / denom
            @. dpolarization_beam[right, :] += bar * leaked * temperature_beam[left, :] / denom
            output_block = @. leaked * factor
            dtemperature_calibration[left] -= sum(bar .* output_block) / temperature_calibration[left]
            dpolarization_calibration[right] -= sum(bar .* output_block) / polarization_calibration[right]
        end

        for (index, (left, right)) in zip(_D1_EE_INDICES, _D1_EE_PAIRS)
            bar = @view EE_bar[left, right, :]
            tt_left, tt_right = minmax(left, right)
            denom = polarization_calibration[left] * polarization_calibration[right]
            factor = @. polarization_beam[left, :] * polarization_beam[right, :] / denom
            TT_block = @view TT[tt_left, tt_right, :]
            TE_left = @view TE[left, right, :]
            TE_right = @view TE[right, left, :]
            EE_block = @view EE[left, right, :]
            leakage = @. gamma[left] * TE_left + gamma[right] * TE_right + gamma[left] * gamma[right] * TT_block
            leaked = @. EE_block + leakage
            @. dEE[left, right, :] += bar * factor
            @. dTE[left, right, :] += bar * factor * gamma[left]
            @. dTE[right, left, :] += bar * factor * gamma[right]
            @. dTT[tt_left, tt_right, :] += bar * factor * gamma[left] * gamma[right]
            @. dgamma[left] += bar * factor * (TE_left + gamma[right] * TT_block)
            @. dgamma[right] += bar * factor * (TE_right + gamma[left] * TT_block)
            @. dpolarization_beam[left, :] += bar * leaked * polarization_beam[right, :] / denom
            @. dpolarization_beam[right, :] += bar * leaked * polarization_beam[left, :] / denom
            output_block = @. leaked * factor
            dpolarization_calibration[left] -= sum(bar .* output_block) / polarization_calibration[left]
            dpolarization_calibration[right] -= sum(bar .* output_block) / polarization_calibration[right]
        end

        parameter_bar = zeros(eltype(output_bar), length(parameters))
        for index in 1:3
            parameter_bar[index] = sum(dgamma[index] .* @view(model.t2p_kernels[index, :]))
        end
        for index in 1:3
            beta = parameters[12 + index]
            main_beam = @view model.main_temperature_beams[index, :]
            normalization = main_beam[model.ell_800_index] + beta * (1 - main_beam[model.ell_800_index])
            numerator = @. main_beam + beta * (temperature_beam[index, :] - main_beam)
            dnumerator = dpolarization_beam[index, :] ./ normalization
            dnormalization = -sum(dpolarization_beam[index, :] .* numerator) / normalization^2
            parameter_bar[12 + index] = sum(dnumerator .* (temperature_beam[index, :] - main_beam)) +
                                        dnormalization * (1 - main_beam[model.ell_800_index])
            @. dtemperature_beam[index, :] += dnumerator * beta
        end
        for mode in 1:9
            parameter_bar[3 + mode] = sum(dtemperature_beam .* @view(model.beam_modes[:, :, mode]))
        end

        dtemperature_calibration .+= dpolarization_calibration .* parameters[19] .* (parameters[20], 1, parameters[21])
        parameter_bar[19] = sum(dpolarization_calibration .* temperature_calibration .* (parameters[20], 1, parameters[21]))
        parameter_bar[20] = dpolarization_calibration[1] * temperature_calibration[1] * parameters[19]
        parameter_bar[21] = dpolarization_calibration[3] * temperature_calibration[3] * parameters[19]
        parameter_bar[16] = dtemperature_calibration[1] * parameters[17] + dtemperature_calibration[2] + dtemperature_calibration[3] * parameters[18]
        parameter_bar[17] = dtemperature_calibration[1] * parameters[16]
        parameter_bar[18] = dtemperature_calibration[3] * parameters[16]

        preinstrument_bar = project_preinstrument(_pack_d1_spectra(dTT, dTE, dEE))
        return NoTangent(), NoTangent(), preinstrument_bar, project_parameters(parameter_bar)
    end

    return output, fast_instrument_pullback
end

function ChainRulesCore.rrule(::typeof(_grids_match), a::AbstractVector, b::AbstractVector)
    return _grids_match(a, b), _ -> (NoTangent(), NoTangent(), NoTangent())
end

end
