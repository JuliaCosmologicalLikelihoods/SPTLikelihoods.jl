using Test
using LinearAlgebra
using NPZ
using SPTLikelihoods
using ADTypes: AutoFiniteDifferences, AutoForwardDiff, AutoMooncake
using DifferentiationInterface
using FiniteDifferences
using ForwardDiff
using Mooncake
using ChainRulesCore
using DelimitedFiles

include("test_helpers.jl")


include("test_spt3g_d1_data.jl")
include("test_spt3g_d1_parameters.jl")
include("test_spt3g_d1_likelihood.jl")
include("test_spt3g_d1_reference.jl")
include("test_spt3g_d1_foregrounds.jl")
include("test_spt3g_d1_instrument.jl")
include("test_spt3g_d1_multipoint.jl")
include("test_spt3g_d1_priors.jl")
include("test_spt3g_d1_ad.jl")
include("test_readme_example.jl")
include("test_benchmark_docs.jl")
include("test_legacy_2018.jl")
include("test_turing.jl")
