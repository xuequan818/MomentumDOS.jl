using TBG_DFT
using LinearAlgebra
using Plots

gauss = [Gaussian(7, 0.05), Gaussian(5, 0.05)]
L = 1
ϵ = 0
model = TbgToy(L, ϵ, gauss)

EcL = 80
EcW = 50
basis = Basis(EcL, EcW, model);

σ = 0.4
xs = collect(-8:0.1:34)
h = 0.1
@time dos = compute_dos_shift_kpm(xs, Gauss(σ), basis, h; Ktrunc = 30);

P = plot(title="ϵ = $ϵ", xs, dos)