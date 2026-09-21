module SPTLikelihoodsMooncakeExt

using SPTLikelihoods: _fixed_covariance_solve, _fast_instrument_core, _grids_match, SPT3GD1ForegroundModel
using LinearAlgebra: Cholesky
using Mooncake: @from_chainrules, MinimalCtx

@from_chainrules MinimalCtx Tuple{typeof(_fixed_covariance_solve), Cholesky{Float64, Matrix{Float64}}, Vector{Float64}}
@from_chainrules MinimalCtx Tuple{typeof(_fast_instrument_core), SPT3GD1ForegroundModel, Vector{Float64}, Vector{Float64}}
@from_chainrules MinimalCtx Tuple{typeof(_grids_match), Vector{Int}, Vector{Int}}

end
