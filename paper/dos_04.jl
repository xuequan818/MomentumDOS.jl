cd("TBG_DFT.jl")
using Pkg
Pkg.activate(".")

#-----------------------------------------------------------------------
using TBG_DFT
using LinearAlgebra
using JLD2

# model setting
gauss = [Gaussian(7, 0.05), Gaussian(5, 0.05)]
L = 1
ϵ = 0 # choose ϵ in {0, 1e-5, -1e-5} to generate DoS for computing the derivatives
model = TbgToy(L, ϵ, gauss)

# initial step
@time dos_int = compute_dos_shift_kpm(collect(1:0.1:2), Gauss(0.5), Basis(20, 10, model), 0.1; M=10, Ktrunc=2);

# basis setting
EcL = 4000
EcW = 80
basis = Basis(EcL, EcW, model);

# DoS computing
σ = 0.4 # Gaussian parameter for test function
h = 0.01
xs = collect(-20:0.1:20) # energy range
@time dos = compute_dos_shift_kpm(xs, Gauss(σ), basis, h; Ktrunc=40, tol=5e-7)

# save data
if iszero(ϵ)
    filename = "dos_0_04.jld2"
elseif ϵ == 1e-5
    filename = "dos_-5_04.jld2"
elseif ϵ == -1e-5
    filename = "dos_-5_04m.jld2"
end
jldsave(filename; ϵ, σ, ER, dos)
