import Pkg
Pkg.activate(@__DIR__)
Pkg.instantiate()

using DelimitedFiles: readdlm
using Plots
include("ratinterptest2.jl")

# Obtain precomputed reflection coefficient data:
filedat = readdlm("s11data.csv", ',', Float64, skipstart = 1)
FGHz = filedat[:, 1]
s11_exact = filedat[:, 2] + im * filedat[:, 3]


knot_indices = vcat([1], length(FGHz) ÷ 5 * (1:5))
errest_max = 1.0

s11_interpolated = zeros(ComplexF64, length(s11_exact))
y_interpolated = copy(s11_interpolated)
s11db_exact = 10*log10.(abs2.(s11_exact))
errest = zeros(length(s11_exact))

stoy(x) = (1-x)/(1+x) # for converting to/from normalized admittance
#stoy(x) = x # for retaining s11
yexact = stoy.(s11_exact)

errest_max_lim_db = -85.0
errest_max_lim = 10^(errest_max_lim_db/20)
anim = @animate while errest_max > errest_max_lim
    global x0j = FGHz[knot_indices]
    global Sj0 = yexact[knot_indices]
    for i in eachindex(y_interpolated, errest, FGHz)
        (y_interpolated[i], errest[i]) = interp_path2(FGHz[i], x0j, Sj0)
    end
    s11db_interpolated = 10*log10.(abs2.(stoy.(y_interpolated)))
    global p = plot(xlim=(FGHz[1],FGHz[end]), ytick=-100:10:0, ylim=(-80,-30), 
    xlabel="Frequency (GHz)", ylabel = "Magnitude (dB)",
    title="Ynorm Interp: $(length(knot_indices)) Knots, $(errest_max_lim_db) dB error criterion", 
    legend=:bottomright)

    plot!(p, FGHz, s11db_exact, color=:blue, label="S11 Exact")
    plot!(p, FGHz, s11db_interpolated, color=:red, label="S11 Interpolated")
    plot!(p, FGHz, 20*log10.(errest), color=:green, label="Error Estimate")
    scatter!(p,  FGHz[knot_indices], fill(-80.0, length(knot_indices)),
    label="Interpolation Knot", markershape = :vline, msize = 3, mcolor = :black)

    global errest_max
    (errest_max, imax) = findmax(errest)
    sort!(push!(knot_indices, imax))
    @show errest_max, length(knot_indices)
end

 
gif(anim, "Ynorm_fitting_animation.gif", fps = 1.5)

scatter!(p,  FGHz[knot_indices], fill(-50.0, length(knot_indices)),
    label=nothing, markershape = :vline, msize = 4, mcolor = :black,
    ylim=(-50,-35), ytick=-50:5:0)
