cd("TBG_DFT.jl")
using Pkg
Pkg.activate(".")

#-----------------------------------------------------------------------

using TBG_DFT
using LinearAlgebra
using JLD2
using FoldsThreads, Folds

gauss = [Gaussian(7, 0.05), Gaussian(5, 0.05)]
L = 1

#initial step
@time dos_int = compute_dos_shift_kpm(collect(1:0.1:2), Gauss(0.5), TbgToy(L, 0, gauss), 20, 10, 0.1; M=10, Ktrunc=2)

# parameters setting
EcL = 6000
EcW = 60
σ = 0.04
h = 0.005
Ktrunc = 24
ER = collect(-16:0.01:18)
incom = 1 + sqrt(2) - round(sqrt(2), digits=5)
ϵs = collect(0.001:0.00016:0.012) .* incom

dos = Vector{Vector{Float64}}(undef, length(ϵs))
Folds.foreach(1:length(ϵs), WorkStealingEx()) do i
	modeli = TbgToy(L, ϵs[i], gauss)
	@time dos[i] = compute_dos_shift_kpm(ER, Gauss(σ), modeli, EcL, EcW, h; Ktrunc,tol=1e-4)
end

jldsave("dos_os_004.jld2"; σ, ϵs, ER, dos)