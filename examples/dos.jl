using TBG_DFT
using LinearAlgebra
using Plots, LaTeXStrings, Plots.Measures

gauss = [Gaussian(7, 0.05), Gaussian(5, 0.05)]
L = 1
EcL = 800
EcW = 40
σ = 0.4
xs = collect(-8:0.1:34);
h = 0.04
ϵ = 0.0
model = TbgToy(L, ϵ, gauss);

basis = Basis(EcL, EcW, model);
@time dos = compute_dos_shift_kpm(xs, Gauss(σ), basis, h; Ktrunc=20, tol=1e-6)

plot(xs, dos, label="ϵ=$ϵ")
plot!(xs, dos0, label="ϵ=0")
plot!(xlims=(-15,-13))

plot!(xlims=(-35,-30))
plot(xs, (dos-dos0)/0.01)
plot(xs,dd)

df = dos - dos0 - 0.01.*dd
plot(xs,df)

using JLD2
f0 = jldopen("dos_0_k_02_1b.jld2", "r")
f1 = jldopen("dos_-6_k_02_1b.jld2", "r")

xs = f0["xs"]
dos0 = f0["dos"]
plot(xs,dos0)

dos1 = f1["dos"]
dd = (dos1 - dos0) ./ (f1["ϵ"] - f0["ϵ"])
plot!(xs, dd)
