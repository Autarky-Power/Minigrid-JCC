# ======================================================================
# Expected.jl
# This code is used to define a minigrid dispatch optimization model 
# with expected values.
# ======================================================================

using JuMP, GLPK, Ipopt, LinearAlgebra, DataFrames
using Plots, XLSX, StatsPlots, Plots.PlotMeasures
using Distributions, Random

sims = 100 # values can be 10 or 100, for getting the right forecasts
Random.seed!(1234)

create_plots = true 

# Set time steps where an outage can occur
T = [i for i=1:24]
# Set outage duration, nut not considered in model
outage = 3
q = 1
prob = 0.9

begin
	dir = @__DIR__
	include("$dir/Utils.jl")
	include("$dir/params.jl")
	season = 1
	inputs_ind.s_max = [mod(i, length(T)) in [0, 1, 2, 3, 4, 5, 6, 22, 23] ? 0.0 : inputs_ind.s_max[i] for i in 1:length(T)+inputs_ind.κ]
end

begin
	# Define the model
	m = Model(Ipopt.Optimizer)
	set_optimizer_attribute(m, "warm_start_init_point", "yes")
	set_optimizer_attribute(m, "linear_solver", "ma27")

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
		@variable(m, E[1:length(T)+inputs.κ] >= 0., start = 0.00)
	end

	# Define objective function
	@objective(m, 
	Min, 
	(sum(inputs_ind.c_Δ .* Δ
			+ inputs_ind.c_b * (b_plus + b_minus)
			- inputs_ind.p_d * inputs_ind.d)
			+ (sum([sum([inputs_ind.c_g_plus[t] * g_plus[t] - inputs_ind.p_g_minus[t] * g_minus[t] + E[t]
									for t in 1:length(T)+inputs_ind.κ if t ∉ τ:τ+inputs_ind.κ])
									for τ in 1:length(T)]) 
				+ (sum([sum([inputs_ind.c_Δ .* r_Δ[t] + inputs_ind.c_b .* r_b_minus[t] 
									for t in τ:τ+inputs_ind.κ]) 
									for τ in eachindex(T)]))) * (inputs_ind.ω  / (length(T)+inputs_ind.κ)) 
			+ sum([inputs_ind.c_g_plus[t] * g_plus[t] - inputs_ind.p_g_minus[t] * g_minus[t] + E[t]
									for t in 1:length(T)+inputs_ind.κ]) * (1-inputs_ind.ω)
		) # * 0.001 # scaling objective function to help solver
		)
	# Define the constraints
	begin
		@constraint(m, 0 .<= Δ + r_Δ .<= inputs_ind.Δ_max)
		@constraint(m, 0 .<= b_minus + r_b_minus .<= inputs_ind.b_max)
		@constraint(m, [t = 2:length(T)+inputs_ind.κ], SOC[t] == SOC[t-1] + (inputs_ind.η_plus * b_plus[t] - inputs_ind.η_minus * b_minus[t]) * Δt / b_cap)
		@constraint(m, SOC[1] == inputs_ind.SOC_0 + (inputs_ind.η_plus * b_plus[1] - inputs_ind.η_minus * b_minus[1]) * Δt / b_cap)
		@constraint(m, c2_vector[i = 1:length(T)+inputs_ind.κ], r_b_minus[i] + r_Δ[i] + inputs.s_max[i] + Δ[i] + b_minus[i] - b_plus[i] - inputs.d[i] >= 0.00)
		@constraint(m, SOC[length(T)] == inputs_ind.SOC_0)
		for t in 1:length(T)+inputs_ind.κ
			for τ in 1:max(t-inputs_ind.κ, 1)
				@constraint(m, inputs_ind.SOC_min <= SOC[t] - sum(inputs_ind.η_minus * r_b_minus[s] * Δt / b_cap for s in τ:min(τ+inputs_ind.κ, t))
							)
			end 
		end
		@constraint(m, [t = 1:length(T)+inputs_ind.κ], 
						SOC[t] .- inputs_ind.SOC_max .<= 0.00
					)
	end
	function define(t, inputs)
		c = inputs_ind.c_g_plus[t] + c_γ_plus
		var = σ[t]
		ψ(y...) = c * var * pdf(Normal(), y[t]/var) + c * y[t] * cdf(Normal(), y[t]/var)
		function ∇ψ(g::AbstractVector{T}, y::T...) where {T}
			for i in eachindex(y)
				if i == t
					g[i] = c * cdf(Normal(), y[i]/var) * 1
				else
					g[i] = 0.00000
				end
			end
			return 			
		end
		return ψ, ∇ψ
	end
	y = inputs_ind.d - inputs_ind.s_max - Δ - g_plus - b_minus + b_plus + g_minus
	for t in eachindex(y)
		register(m, Symbol("expected_$t"), length(y), define(t, inputs)[1], define(t, inputs)[2])
		add_nonlinear_constraint(m, :($(Symbol("expected_$t"))($(y...)) == $(E[t])))
	end
end

# Solve
@time JuMP.optimize!(m) 
solution_summary(m, verbose=true)

# If optimal, save results to excel
if termination_status(m) in [MOI.LOCALLY_SOLVED, MOI.OPTIMAL]
	# Write results to excel file
	dump_results_expected(T, m, season, inputs_ind)
else
	println(termination_status(m))
end

if create_plots
	# Plot of dispatch
	begin
		if termination_status(m) in [MOI.LOCALLY_SOLVED, MOI.OPTIMAL] 
			gr(size=(800,600))
			direct_solar = [inputs_ind.s_max[i] == 0.00 ? 0.00 : inputs_ind.d[i]-value.(m[:g_plus])[i]-value.(m[:b_minus])[i] for i in 1:length(T)+inputs_ind.κ]
			ticklabel = [1:length(T)+inputs_ind.κ]
			d = groupedbar([value.(m[:Δ]) direct_solar value.(m[:b_minus]) value.(m[:g_plus]) -value.(m[:g_minus]) -value.(m[:b_plus])],
					bar_position = :stack,
					bar_width=0.7,
					xticks=(1:2:length(T)+inputs_ind.κ),
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
			# title!("\n\nopt dispatch for the expected-value model, p = $(inputs_ind.p), objective value = $(-round(objective_value(m), digits=3)) €", titlelocation=:center) #, fontsize=3)
			plot!(inputs_ind.d, label = "forecast. demand", linewidth=2, color=:black)
			plot!(inputs_ind.s_max, label = "forecast. total solar", linewidth=2, color=:orange)
		end 
	end 
	savefig(d, "ExpectedModel_$(inputs_ind.p *100).png")

	# Plot of reserves
	begin
		if termination_status(m) in [MOI.LOCALLY_SOLVED, MOI.OPTIMAL] 
			gr(size=(800,600))
			ticklabel = [1:length(T)+inputs_ind.κ]
			res = groupedbar([value.(m[:r_b_minus]) value.(m[:r_Δ])],
					bar_position = :stack,
					bar_width=0.7,
					xticks=(1:2:length(T)+inputs_ind.κ),
					tickfontsize = 11,
					ylims=(0, 20),
					label=["battery discharge reserve" "diesel reserve"],
					legend=:topleft,
					legendfontsize=12,
					top_margin = 10mm,
					xlabel = "\n\nTime in h",
					ylabel = "\n\nPower in kW",
					labelfontsize = 13,
					color = [:magenta4 :thistle])
			# title!("\n\nopt reserves for the expected-value model, p = $(inputs_ind.p), objective value = $(-round(objective_value(m), digits=3)) €", titlelocation=:center) #, fontsize=3)
		end 
	end
	savefig(res, "Reserves_ExpectedModel_$(inputs_ind.p *100).png")

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
			# title!("\n\nopt battery dispatch for the expected-value model, p = $(inputs_ind.p), objective value = $(-round(objective_value(m),digits=3)) €", titlelocation=:center) #, fontsize=3)
			plot!(twinx(), value.(m[:SOC]), ylims=(15, 80), label = "SOC", linewidth=2, color=:black, ylabel="\n\nSOC in %\n\n", legendfontsize=12, tickfontsize = 14)
		end 
	end
	savefig(batt, "Battery_ExpectedModel_$(inputs_ind.p *100).png")
end