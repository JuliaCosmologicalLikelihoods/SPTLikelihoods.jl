using DelimitedFiles
using SPTLikelihoods

const SPT3G_D1_FIXTURE_DIR = joinpath(@__DIR__, "fixtures", "spt3g_d1_full")
const SPT3G_D1_TEST_LIKE_DIR = joinpath(@__DIR__, "fixtures", "spt3g_d1_test_like")

_read_matrix(path::AbstractString) = readdlm(path, Float64; comments=true, comment_char='#')

function _load_txt_parameters(path::AbstractString)
    raw = readdlm(path, String; comments=true, comment_char='#')
    return Dict{String,Float64}(raw[i, 1] => parse(Float64, raw[i, 2]) for i in 1:size(raw, 1))
end

function _load_test_models(fixture_dir::AbstractString=SPT3G_D1_FIXTURE_DIR)
    params_dict = _load_txt_parameters(joinpath(fixture_dir, "baseline_parameters.txt"))
    params = SPT3GD1Parameters(params_dict)
    cmb_mat = _read_matrix(joinpath(fixture_dir, "cmb_spectra.txt"))
    ells = Int.(cmb_mat[:, 1])
    cmb = SPT3GD1CMBTheory(ells, cmb_mat[:, 2], cmb_mat[:, 3], cmb_mat[:, 4])
    model = SPT3GD1ForegroundModel(fixture_dir)
    return (params=params, params_dict=params_dict, cmb=cmb, model=model)
end
