cd("MomentumDOS.jl")
using Pkg
Pkg.activate(".")

using MomentumDOS
using LinearAlgebra
using JLD2

function test_dos(ϵ::Float64, σ::Float64)
	if σ == 0.4
		EcL = 4000
		EcW = 80
		h = 0.05
	elseif σ == 0.08
		EcL = 5000
		EcW = 80
		h = 0.01
	elseif σ == 0.04
		EcL = 6000
		EcW = 80
		h = 0.005
	end
	
    gauss = [Gaussian(7, 0.05), Gaussian(5, 0.05)]
    L = 1
    model = TbgToy(L, ϵ, gauss)

    @time dos_int = compute_dos_shift_kpm(collect(1:0.1:2), Gauss(0.5), Basis(20, 10, model), 0.1; M=10, Ktrunc=2)

    basis = Basis(EcL, EcW, model)
    ER = collect(-35:0.1:35)
    @time dos = compute_dos_shift_kpm(ER, Gauss(σ), basis, h; Ktrunc=40)

	ϵsr = string(ϵ)[3:end]
    if occursin('e', ϵsr)
		ϵsr = "-5"
	end

    σsr = replace(string(σ), "." => "")

	if ϵ >= 0
		file_name = "dos_$(ϵsr)_$(σsr).jld2"
	else
		file_name = "dos_$(ϵsr)_$(σsr)m.jld2"
	end

    jldsave(file_name; ϵ, σ, xs=ER, dos)
end
