 using JLD2, Plots, Plots.Measures
 using LinearAlgebra

ϵ = -4
ref = jldopen("test_005_sigma/ref_$(ϵ)_005.jld2", "r")
#Wtest = jldopen("test_results_LW/Wtest_$(ϵ)_01.jld2", "r")
Ltest = jldopen("test_005_sigma/Ltest_$(ϵ)_005.jld2", "r")

#e2W = [norm(x - ref["ldos_ref"]) for x in Wtest["ldosTest"]]
#eInfW = [norm(x - ref["ldos_ref"],Inf) for x in Wtest["ldosTest"]]

e2L = [norm(x - ref["ldos_ref"]) for x in Ltest["ldosTest"]]
eInfL = [norm(x - ref["ldos_ref"][:,1:21], Inf) for x in Ltest["ldosTest"]]


P1 = plot(Wtest["WTest"], e2W, yscale=:log10, ylabel="", xlabel="W", guidefontsize=22, color=:black, title="ϵ = $(ϵ)", label="2", tickfontsize=20, legendfontsize=20, legend=:topright, grid=:off, box=:on, size=(740, 600), titlefontsize=30, left_margin=2mm, right_margin=2mm, top_margin=3mm, lw=2, marker=:circle, markersize=8, markercolor=:white, markerstrokecolor=:black)
plot!(P1, Wtest["WTest"], eInfW, label="∞", color=:red, lw=2, marker=:circle, markersize=8, markercolor=:white, markerstrokecolor=:red)


P2 = plot(Ltest["LTest"], e2L, yscale=:log10, ylabel="", xlabel="L", guidefontsize=22, color=:black, title="", label="2", tickfontsize=20, legendfontsize=20, legend=:topright, grid=:off, box=:on, size=(740, 600), titlefontsize=30, left_margin=2mm, right_margin=4mm, top_margin=3mm, lw=2, marker=:circle, markersize=8, markercolor=:white, markerstrokecolor=:black)
plot!(P2, Ltest["LTest"], eInfL, label="∞", color=:red, lw=2, marker=:circle, markersize=8, markercolor=:white, markerstrokecolor=:red)

P = plot([P1, P2]..., size=(740, 900), ylabel="Error",left_margin=6mm,layout = grid(2, 1, heights=[0.5, 0.5]))
#savefig("error_lods_-4.png")
#"W = 120, L = 1600 for ϵ = 0"