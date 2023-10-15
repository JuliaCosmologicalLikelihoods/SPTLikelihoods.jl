# SPTLikelihoods.jl

Repository containing a Julia native likelihood of the SPT Likelihoods.

In this moment we support:

- SPT-3G 2018 TT/TE/EE, as presented in [Balkenhol et al. (2022)](https://arxiv.org/abs/2212.05642)

We plan to include more likelihoods in the forthcoming future.

Our likelihoods, if coupled with a differentiable Boltzmann solver or emulator, can be used with gradient based samplers such as HMC, NUTS, or MCHMC.

## Usage

 Usage of the code in this repository is warranted, upon citing the relevant papers for the data employed and the paper when these likelihoods were presented:

 - Bonici, Bianchini, and Ruiz-Zapatero, _Capse.jl: efficient and auto-differentiable CMB power spectra emulation_ [![arXiv](https://img.shields.io/badge/arXiv-2307.14339-b31b1b.svg)](https://arxiv.org/abs/2307.14339)
