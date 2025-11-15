# TBG_DFT.jl

A Julia packge for computing the density of states (DoS) of the 1D toy [incommensurate Hamiltonian](https://arxiv.org/pdf/2510.15369) using the [momentum-space](https://epubs.siam.org/doi/abs/10.1137/23M1553650) method.

## Installation
TBG_DFT.jl is an unregistered package and therefore needs to be downloaded or cloned to the user's local computer first, and then installed by running

```julia
julia> cd("your-local-path/TBG_DFT.jl")
julia> using Pkg
julia> Pkg.activate(".")
julia> Pkg.instantiate()
```

## Usage
Here are some examples showing how to use the packge.

Simulation of a 1D toy model:  
```math
	H_{\epsilon} = -\frac{1}{2}  \frac{d^2}{dx^2} + V(x,(1+\epsilon)x)
```
where the potentials are given by
```math
\begin{align*}
&V(x,(1+\epsilon)x) := \sum_{R\in \mathbb{Z}} \Big(v_1\big(x-R\big)+v_2\big((1+\epsilon)x-R\big)\Big),\\
& v_j(x)=-\frac{A_j}{\sqrt{2\pi\sigma_j^2}} e^{-\frac{|x|^2}{2\sigma_j^2}},\qquad j= 1,2.
\end{align*}
```
```julia
using TBG_DFT
using Plots

# Define the 1D model
gauss = [Gaussian(7, 0.05), Gaussian(5, 0.05)]
L = 1
ϵ = pi/200
model = TbgToy(L, ϵ, gauss);

# Define the basis
EcL = 300
EcW = 20
basis = Basis(EcL, EcW, model);

# Compute the DoS
ER = collect(-8:0.1:34) # Energy range
h = 0.1 
σ = 0.4 # Gaussian parameter for test function
@time dos = compute_dos_shift_kpm(ER, Gauss(σ), basis, h; Ktrunc=20, tol=1e-4);

plot(ER, dos, label="ϵ=$ϵ")
```
