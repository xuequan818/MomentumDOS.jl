# TBG_DFT.jl

A Julia packge for computing the density of states (DoS) of the [1D toy incommensurate Hamiltonian](https://arxiv.org/pdf/2510.15369) using the [momentum space](https://epubs.siam.org/doi/abs/10.1137/23M1553650) method.

## Installation
TBG_DFT.jl is an unregistered package and therefore needs to be downloaded or cloned to the user's local computer first, and then installed by running

```julia
julia> cd("your-local-path/GeneralizedBM.jl")
julia> using Pkg
julia> Pkg.activate(".")
julia> Pkg.instantiate()
```
