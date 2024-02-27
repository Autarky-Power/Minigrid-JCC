# ======================================================================
# RegularProb.jl
# This code is used to define a minigrid dispatch optimization model 
# with the regular "deterministic" model.
# ======================================================================
using DataFrames
using Distributions
using GLPK
using JuMP
using Random
using LinearAlgebra
using Plots, XLSX, StatsPlots, Plots.PlotMeasures
using XLSX

# Fix the seed
Random.seed!(1234)

create_plots = true

begin
    # Set time steps where an outage can occur
    T = [i for i=1:24]
    # Set outage duration
    outage = 3
    # Scale demand 
    q = 1.
    # Determine season: values can be 1 for "dry" or 2 for "rainy"
    season = 1 
    # Set number of simulations for getting the right forecasts. Values can be 10 or 100. 
    sims = 100
	 
	prob = 0.9
end 

begin
	dir = @__DIR__
	include("$dir/Utils.jl")
end
include("$dir/params.jl")
inputs_ind.s_max = [mod(i, length(T)) in [0, 1, 2, 3, 4, 5, 6, 22, 23] ? 0.0 : inputs_ind.s_max[i] for i in 1:length(T)+inputs_ind.κ]


begin
	# Define the model
	m = Model(GLPK.Optimizer)

	# Define positive decision variables
	@variable(m, 0 <= Δ[1:length(T)+inputs_ind.κ] <= inputs_ind.Δ_max)
	@variable(m, 0 <= b_plus[1:length(T)+inputs_ind.κ] <= inputs_ind.b_max)
	@variable(m, 0 <= b_minus[1:length(T)+inputs_ind.κ] <= inputs_ind.b_max)
	@variable(m, 0 <= g_plus[i = 1:length(T)+inputs_ind.κ] <= inputs_ind.g_max[i])
	@variable(m, 0 <= g_minus[i = 1:length(T)+inputs_ind.κ] <= inputs_ind.g_max[i])
	@variable(m, inputs_ind.SOC_min <= SOC[1:length(T)+inputs_ind.κ] <= inputs_ind.SOC_max, start = inputs_ind.SOC_0)
	# @variable(m, 0 <= s[i = 1:length(T)+inputs_ind.κ] <= inputs_ind.s_max[i])

	@objective(m, Min, sum(inputs_ind.c_Δ .* Δ +
		                   inputs_ind.c_g_plus .* g_plus - 
						   inputs_ind.p_g_minus .* g_minus -
		  			       inputs_ind.p_d * inputs_ind.d +
		  			       inputs_ind.c_b * (b_plus + b_minus)
						   )
		            )

	# Define the constraints
	@constraint(m, [t = 2:length(T)+inputs_ind.κ], SOC[t] == SOC[t-1] + (inputs_ind.η_plus * b_plus[t] - inputs_ind.η_minus * b_minus[t]) * Δt / b_cap)
	@constraint(m, SOC[1] == inputs_ind.SOC_0 + (inputs_ind.η_plus * b_plus[1] - inputs_ind.η_minus * b_minus[1]) * Δt / b_cap)
	@constraint(m, c1_vector, inputs_ind.s_max + Δ + g_plus + b_minus - b_plus - g_minus - inputs_ind.d .== 0)
	@constraint(m, SOC[length(T)] == inputs_ind.SOC_0)
end

# Solve
@time JuMP.optimize!(m)
solution_summary(m, verbose=true)

if termination_status(m) in [MOI.LOCALLY_SOLVED, MOI.OPTIMAL]
	dump_results(T, m, season, inputs_ind)
else
	println(termination_status(m))
end

if create_plots
	begin
		gr(size=(800,600))
		direct_solar = [inputs_ind.s_max[i] == 0.00 ? 0.00 : inputs_ind.s_max[i]-value.(m[:g_minus])[i]-value.(m[:b_plus])[i] for i in 1:length(T)+inputs_ind.κ]
		ticklabel = [1:length(T)+inputs_ind.κ]
		d = groupedbar([value.(m[:Δ]) direct_solar value.(m[:b_minus]) value.(m[:g_plus]) -value.(m[:g_minus]) -value.(m[:b_plus])],
				bar_position = :stack,
				bar_width=0.7,
				xticks=(1:2:length(T)+inputs_ind.κ), #, ticklabel),
				ylims=(-20,21),
				label=["diesel" "direct solar" "battery discharge" "grid import" "grid export" "battery charge"],
				legend=:bottomright,
				top_margin = 10mm,
				xlabel = "\n\nTime in h",
				ylabel = "\n\nPower in kW",
				labelfontsize = 13,
				legendfontsize=10,
				tickfontsize = 11,
				color= [:grey :orange :green4 :blue4 :royalblue1 :chartreuse2])
		# title!("\n\nopt dispatch for regular model, objective value = $(-round(objective_value(m), digits=3)) €", titlelocation=:center) #, fontsize=3)
		plot!(inputs_ind.d, label = "forecast. demand", linewidth=2, color=:black)
		plot!(inputs_ind.s_max, label = "forecast. total solar", linewidth=2, color=:orange)
	end 
	savefig("RegularModel.png")

	# Plot of battery profile
	begin
		if termination_status(m) in [MOI.LOCALLY_SOLVED, MOI.OPTIMAL] 
			gr(size=(900,700))
			ticklabel = [1:length(T)+inputs_ind.κ]
			batt = groupedbar([value.(m[:b_minus]) -value.(m[:b_plus])],
					bar_position = :stack,
					bar_width=0.7,
					xticks=(1:2:length(T)+inputs_ind.κ),
					ylims=(-20, 20),
					label=["battery discharge" "battery charge"],
					legend=:topleft,
					top_margin = 10mm,
					xlabel = "\n\nTime in h",
					ylabel = "\n\nPower in kW",
					labelfontsize = 13,
					legendfontsize=12,
					tickfontsize = 14,
					color = [:green4 :chartreuse2])
			# title!("\n\nopt battery dispatch for regular model, objective value = $(-round(objective_value(m), digits=3)) €", titlelocation=:center) #, fontsize=3)
			plot!(twinx(), value.(m[:SOC]), ylims=(15, 80), label = "SOC", linewidth=2, color=:black, ylabel="\n\nSOC in %\n\n", legendfontsize=12, tickfontsize = 14)
		end 
	end
	savefig("Battery_RegularModel.png")
end