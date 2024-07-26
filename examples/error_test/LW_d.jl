 using JLD2, Plots, Plots.Measures
 using LinearAlgebra

ldos0 = jldopen("test_results_LW/Wtest_0_01.jld2", "r")
ldos1 = jldopen("test_results_LW/Wtest_-4_01.jld2", "r")

dldos1 = (ldos1["ldosTest"][1] .- ldos0["ldosTest"][1]) ./ (2pi * (ldos1["ϵ"] - ldos0["ϵ"]))
dldos2 = (ldos1["ldosTest"][end] .- ldos0["ldosTest"][end]) ./ (2pi * (ldos1["ϵ"] - ldos0["ϵ"]))
@show norm(dldos1 - dldos2)
heatmap(ldos1["xs"], ldos1["ys"], ldos1["ldosTest"][end])
