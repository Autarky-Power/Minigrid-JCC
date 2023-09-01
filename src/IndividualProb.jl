# ======================================================================
# IndividualProb.jl
# This code is used to define a minigrid dispatch optimization model 
# with Individual Probabilistic (Chance) Constraints.
# ======================================================================

using JuMP, GLPK, LinearAlgebra, DataFrames
using Plots, XLSX, StatsPlots, Plots.PlotMeasures
using Distributions, Random

sims = 100 # values can be 10 or 100, for getting the right forecasts
Random.seed!(1234)

# Set time steps where an outage can occur
T = [i for i=1:24]
# Set outage duration, nut not considered in model
outage = 3
q = 1
prob = 0.92

begin
	dir = @__DIR__
	include("$dir/Utils.jl")
	include("$dir/params.jl")
	season = 1
	inputs_ind.s_max = [mod(i, length(T)) in [0, 1, 2, 3, 4, 5, 6, 22, 23] ? 0.0 : inputs.s_max[i] for i in 1:length(T)+inputs.κ]
end

begin
	# Define the model
	m = Model(GLPK.Optimizer)

	# Define variables
	begin
		# Define positive decision variables
		@variable(m, 0 <= Δ[1:length(T)+inputs_ind.κ] <= inputs_ind.Δ_max)
		@variable(m, 0 <= b_plus[1:length(T)+inputs_ind.κ] <= inputs_ind.b_max)
		@variable(m, 0 <= b_minus[1:length(T)+inputs_ind.κ] <= inputs_ind.b_max)
		@variable(m, 0 <= g_plus[i = 1:length(T)+inputs_ind.κ] <= inputs_ind.g_max[i])
		@variable(m, 0 <= g_minus[i = 1:length(T)+inputs_ind.κ] <= inputs_ind.g_max[i])
		@variable(m, r_Δ[1:length(T)+inputs_ind.κ] >= 0)
		@variable(m, r_b_plus[1:length(T)+inputs_ind.κ] >= 0)
		@variable(m, r_b_minus[1:length(T)+inputs_ind.κ] >= 0)
		@variable(m, SOC[1:length(T)+inputs_ind.κ] >= inputs_ind.SOC_min, start = inputs_ind.SOC_0)
		@variable(m, 0 <= s[i = 1:length(T)+inputs_ind.κ] <= inputs_ind.s_max[i])
		# @variable(m, ω[1:length(T)+inputs_ind.κ] >= 0)
	end

	# Define objective function
	@objective(m, Min, sum(inputs_ind.c_Δ .* Δ -
		  			       inputs_ind.p_d * inputs_ind.d +
		  			       inputs_ind.c_b * (b_plus + b_minus)) +
						   sum([inputs_ind.c_g_plus[t] * g_plus[t] - inputs_ind.p_g_minus[t] * g_minus[t] for t in 1:length(T)+inputs_ind.κ]) * (1-inputs_ind.ω) +
						   (length(T) == 1 ? 0 : sum([sum([inputs_ind.c_g_plus[t] * g_plus[t] - inputs_ind.p_g_minus[t] * g_minus[t] for t in 1:length(T)+inputs_ind.κ if t ∉ τ]) + 
						   inputs_ind.c_Δ .* r_Δ[τ] + inputs_ind.c_b * r_b_minus[τ] for τ in 1:length(T)+inputs_ind.κ]) * (inputs_ind.ω  / (length(T)+inputs_ind.κ)) + 
						   sum([inputs_ind.c_g_plus[t] * g_plus[t] - inputs_ind.p_g_minus[t] * g_minus[t] for t in 1:length(T)+inputs_ind.κ]) * (1-inputs_ind.ω))
		            )
	# Define the constraints
	begin
		@constraint(m, 0 .<= Δ + r_Δ .<= inputs_ind.Δ_max)
		# @constraint(m, 0 .<= b_plus + r_b_plus .<= inputs_ind.b_max)
		@constraint(m, 0 .<= b_minus + r_b_minus .<= inputs_ind.b_max)
		@constraint(m, [t = 2:length(T)+inputs.κ], SOC[t] == SOC[t-1] + (inputs_ind.η_plus * b_plus[t] - inputs_ind.η_minus * b_minus[t]) * Δt / b_cap)
		@constraint(m, SOC[1] == inputs_ind.SOC_0 + (inputs_ind.η_plus * b_plus[1] - inputs_ind.η_minus * b_minus[1]) * Δt / b_cap)
		# @constraint(m, [t = 2:length(T)+inputs.κ], inputs_ind.SOC_min .<= SOC[t-1] + (inputs_ind.η_plus * (b_plus[t] + r_b_plus[t]) + inputs_ind.η_minus * (b_minus[t] + r_b_minus[t])) * Δt / b_cap .<= inputs_ind.SOC_max)
		# @constraint(m, inputs_ind.SOC_min <= inputs_ind.SOC_0 + (inputs_ind.η_plus * (b_plus[1] + r_b_plus[1]) +
		#                inputs_ind.η_minus * (b_minus[1] + r_b_minus[1])) * Δt / b_cap <= inputs_ind.SOC_max)
		@constraint(m, [t = 2:length(T)+inputs.κ], inputs_ind.SOC_min .<= SOC[t-1] + (inputs_ind.η_plus * (b_plus[t]) + inputs_ind.η_minus * (b_minus[t] + r_b_minus[t])) * Δt / b_cap .<= inputs_ind.SOC_max)
		@constraint(m, inputs_ind.SOC_min <= inputs_ind.SOC_0 + (inputs_ind.η_plus * (b_plus[1]) +
									  inputs_ind.η_minus * (b_minus[1] + r_b_minus[1])) * Δt / b_cap <= inputs_ind.SOC_max)
		@constraint(m, c1_vector, s + Δ + g_plus + b_minus - b_plus - g_minus - inputs_ind.d .== 0)
		# @constraint(m, s + Δ + g_plus + b_minus - b_plus - g_minus - inputs_ind.d - ω .== 0)
		# @constraint(m, [i = 1:length(T)+inputs.κ], r_b_minus[i] - r_b_plus[i] + r_Δ[i] >= quantile(ξ[i], inputs_ind.p))
		@constraint(m, c2_vector[i = 1:length(T)+inputs.κ], r_b_minus[i] + r_Δ[i] - g_plus[i] + g_minus[i] >= quantile(ξ[i], inputs_ind.p))
		@constraint(m, SOC[length(T)] == inputs_ind.SOC_0)
	end
end

# Solve
@time JuMP.optimize!(m) 
solution_summary(m, verbose=true)

# Set number of simulations used for sampling
nsims = 10 
# If optimal, save results to excel
if termination_status(m) in [MOI.LOCALLY_SOLVED, MOI.OPTIMAL]
	# Write results to excel file
	dump_results_ind(T, m, season, inputs_ind)
	# Create the sample
	γ = create_sample(ξ, inputs_ind.κ, nsims)
	# Create and save plots
	savefig(plot_results_individ(T, m, γ, inputs_ind.κ, nsims)[1], "Errors_withR_individualCC_season$(lpad(season, 1))_prob$(lpad(inputs_ind.p, 3))_$(nsims).png")
	savefig(plot_results_individ(T, m, γ, inputs_ind.κ, nsims)[2], "Errors_withoutR_individualCC__season$(lpad(season, 1))_prob$(lpad(inputs_ind.p, 3))_$(nsims).png")
else
	println(termination_status(m))
end

# Plot of dispatch
begin
	if termination_status(m) in [MOI.LOCALLY_SOLVED, MOI.OPTIMAL] 
		gr(size=(800,600))
		ticklabel = [1:length(T)+inputs_ind.κ]
		d = groupedbar([value.(m[:Δ]) value.(m[:s])-value.(m[:g_minus])-value.(m[:b_plus]) value.(m[:b_minus]) value.(m[:g_plus]) -value.(m[:g_minus]) -value.(m[:b_plus])],
				bar_position = :stack,
				bar_width=0.7,
				xticks=(1:1:length(T)+inputs_ind.κ),
				ylims=(-20,21),
				label=["diesel" "direct solar" "battery discharge" "grid import" "grid export" "battery charge"],
				legend=:bottomright,
				top_margin = 10mm,
				xlabel = "\n\nTime in h",
				ylabel = "\n\nPower in kW",
				color= [:grey :orange :green4 :blue4 :royalblue1 :chartreuse2])
		# title!("\n\nopt dispatch for ICC model, p = $(inputs_ind.p), objective value = $(-round(objective_value(m), digits=3)) €", titlelocation=:center) #, fontsize=3)
		plot!(inputs_ind.d, label = "demand", linewidth=2, color=:black)
		plot!(inputs_ind.s_max, label = "total solar", linewidth=2, color=:orange)
	end 
end 
savefig(d, "ICCModel_$(inputs_ind.p *100).png")

# Plot of reserves
begin
	if termination_status(m) in [MOI.LOCALLY_SOLVED, MOI.OPTIMAL] 
		gr(size=(800,600))
		ticklabel = [1:length(T)+inputs_ind.κ]
		res = groupedbar([value.(m[:r_b_minus]) value.(m[:r_Δ])],
				bar_position = :stack,
				bar_width=0.7,
				xticks=(1:1:length(T)+inputs_ind.κ),
				ylims=(0, 20),
				label=["battery discharge reserve" "diesel reserve"],
				legend=:topleft,
				top_margin = 10mm,
				xlabel = "\n\nTime in h",
				ylabel = "\n\nPower in kW",
				color = [:magenta4 :thistle])
		# title!("\n\nopt reserves for ICC model, p = $(inputs_ind.p), objective value = $(-round(objective_value(m), digits=3)) €", titlelocation=:center) #, fontsize=3)
	end 
end
savefig(res, "Reserves_ICCModel_$(inputs_ind.p *100).png")

# Plot of battery profile
begin
	if termination_status(m) in [MOI.LOCALLY_SOLVED, MOI.OPTIMAL] 
		gr(size=(900,700))
		ticklabel = [1:length(T)+inputs_ind.κ]
		batt = groupedbar([value.(m[:b_minus]) -value.(m[:b_plus])],
				bar_position = :stack,
				bar_width=0.7,
				xticks=(1:1:length(T)+inputs_ind.κ),
				ylims=(-20, 20),
				label=["battery discharge" "battery charge"],
				legend=:topleft,
				top_margin = 10mm,
				xlabel = "\n\nTime in h",
				ylabel = "\n\nPower in kW",
				color = [:green4 :chartreuse2])
		# title!("\n\nopt battery dispatch for ICC model, p = $(inputs_ind.p), objective value = $(-round(objective_value(m),digits=3)) €", titlelocation=:center) #, fontsize=3)
		plot!(twinx(), value.(m[:SOC]), label = "state of charge", linewidth=2, color=:black, ylabel="\n\nSOC in %\n\n")
	end 
end
savefig(batt, "Battery_ICCModel_$(inputs_ind.p *100).png")

# # Get shadow prices 
# sp_balancing = [shadow_price(c1_vector[i]) for i in 1:length(T)+inputs_ind.κ]
# sp_CC = [shadow_price(c2_vector[i]) for i in 1:length(T)+inputs_ind.κ]

# println(sp_balancing)
# println(sp_CC)