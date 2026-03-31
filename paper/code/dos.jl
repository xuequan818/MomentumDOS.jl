cd("TBG_DFT.jl")
using Pkg
Pkg.activate(".")

#-----------------------------------------------------------------------
#initial step

using TBG_DFT
using LinearAlgebra
using JLD2


function test_dos(ϵ::Float64; σ="04")
	if σ == "04"
		EcL = 4000
		EcW = 80
		h = 0.05
	elseif σ == "008"
		EcL = 5000
		EcW = 80
		h = 0.01
	elseif σ == "004"
		EcL = 6000
		EcW = 80
		h = 0.005
	end
	
    gauss = [Gaussian(7, 0.05), Gaussian(5, 0.05)]
    L = 1
    model = TbgToy(L, ϵ, gauss)

    @time dos_int = compute_dos_shift_kpm(collect(1:0.1:2), Gauss(0.5), Basis(20, 10, model), 0.1; M=10, Ktrunc=2)

    basis = Basis(EcL, EcW, model)
    xs = collect(-35:0.1:35)
    @time dos = compute_dos_shift_kpm(ER, Gauss(σ), basis, h; Ktrunc=40)

	ϵtr = string(ϵ)[3:end]
    if occursin('e', ϵtr)
		ϵtr = "-5"
	end

	if ϵ >= 0
		file_name = "dos_$(ϵtr)_$(σ).jld2"
	else
		file_name = "dos_$(ϵtr)_$(σ)m.jld2"
	end

    jldsave(file_name; ϵ, σ, xs, dos)
end
