cd("TBG_DFT.jl")
using Pkg
Pkg.activate(".")

#-----------------------------------------------------------------------
using TBG_DFT
using LinearAlgebra
using JLD2

gauss = [Gaussian(7, 0.05), Gaussian(5, 0.05)]
L = 1
incom = 1 + sqrt(2) - round(sqrt(2), digits=5)
ϵ = 0.01 * incom
model = TbgToy(L, ϵ, gauss)

@time dos_int = compute_dos_shift_kpm(collect(1:0.1:2), Gauss(0.5), Basis(20, 10, model), 0.1; M=10, Ktrunc=2);

EcL = 6000
EcW = 80
basis = Basis(EcL, EcW, model);

σ = 0.04
h = 0.005
xs = collect(-20:0.01:20)
@time dos = compute_dos_shift_kpm(xs, Gauss(σ), basis, h; Ktrunc=35, tol=1e-4)

jldsave("dos_01_004.jld2"; ϵ, σ, xs, dos)