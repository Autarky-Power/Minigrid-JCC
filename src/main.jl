# ======================================================================
# main.jl
# This code is used to define apply a scenario to a minigrid dispatch 
# optimization model with Joint Probabilistic (Chance) Constraints.

# It can apply the model with the Genz method or with SRD.
# The file also creates plots to visualize the results.
# ======================================================================

using DataFrames, Distributions, LinearAlgebra, MvNormalCDF, Random, SpecialFunctions
using JuMP, Ipopt, NLopt
using Plots, XLSX, StatsPlots, Plots.PlotMeasures

# Fix the seed
Random.seed!(1234)

create_plots = true

# Import input files, functions and parameters
begin
	dir = @__DIR__
	include("$dir/Utils.jl")
    include("$dir/JointProb.jl")
end

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
	# Set probability level 
	prob = 0.9
end 

# Define scenario
scenario = Scenario(T, outage, q, season)

# Define parameters
include("$dir/params.jl")
inputs.s_max = [mod(i, length(T)) in [0, 1, 2, 3, 4, 5, 6, 22, 23] ? 0.0 : inputs.s_max[i] for i in 1:length(T)+inputs.κ]

# Specify solver
solver = "Ipopt"

# Define model and model options

if solver == "Ipopt"
	m = Model(Ipopt.Optimizer)
	set_optimizer_attribute(m, "warm_start_init_point", "yes")
	set_optimizer_attribute(m, "linear_solver", "ma27")
	set_optimizer_attribute(m, "max_iter", 6000)
elseif solver == "Pavito"
	m = Model(
		optimizer_with_attributes(
		Pavito.Optimizer,
		"mip_solver" => optimizer_with_attributes(GLPK.Optimizer),
		"cont_solver" =>
			optimizer_with_attributes(Ipopt.Optimizer, "print_level" => 0),
			))
elseif solver == "NLopt"
	m = Model(NLopt.Optimizer)
	set_optimizer_attribute(m, "algorithm", :LD_SLSQP)
	set_optimizer_attribute(m, "local_optimizer", :LD_LBFGS)

end 

# Build the model
m = define_model_JCC_Genz(scenario, inputs, m)
# m = define_model_JCC_SRD(scenario, inputs, m)

# Solve
@time JuMP.optimize!(m)
solution_summary(m, verbose=true)

# Extract optimal values of variables
vars = Dict(
    k => value.(v) for 
    (k, v) in object_dictionary(m) if v isa AbstractArray{VariableRef}
)

# Set number of simulations used for sampling
nsims = 10
# Optimality as condition to create output excels & figs
if termination_status(m) in [MOI.LOCALLY_SOLVED, MOI.OPTIMAL]
	dump_results_joint(T, m, season, inputs)
	# Create the sample
	γ = create_sample_joint(ξ, Σ, inputs.κ, nsims)
	# Create and save plots
	for f in 1:length(T)
	    savefig(plot_results(T, m, γ, f, inputs.κ, nsims)[1], "Errors_withR_season$(lpad(season, 1))_kappa$(lpad(inputs.κ, 2, "0"))_prob$(lpad(inputs.p, 3))_window$(lpad(f, 2, "0"))_$(nsims).png")
	    savefig(plot_results(T, m, γ, f, inputs.κ, nsims)[2], "Errors_withoutR__season$(lpad(season, 1))_kappa$(lpad(inputs.κ, 2, "0"))_prob$(lpad(inputs.p, 3))_window$(lpad(f, 2, "0"))_$(nsims).png")
	end
else
	println(termination_status(m))
end

if create_plots
	# Plot of dispatch
	begin
		if termination_status(m) in [MOI.LOCALLY_SOLVED, MOI.OPTIMAL] 
			direct_solar = [inputs.s_max[i] == 0.00 ? 0.00 : inputs_ind.d[i]-value.(m[:g_plus])[i]-value.(m[:b_minus])[i] for i in 1:length(T)+inputs.κ] # inputs.s_max[i]-value.(m[:g_minus])[i]-value.(m[:b_plus])[i] for i in 1:length(T)+inputs.κ]
			gr(size=(800,600))
			ticklabel = [1:length(T)+inputs.κ]
			d = groupedbar([value.(m[:Δ]) direct_solar value.(m[:b_minus]) value.(m[:g_plus]) -value.(m[:g_minus]) -value.(m[:b_plus])],
					bar_position = :stack,
					bar_width=0.7,
					xticks=(1:2:length(T)+inputs.κ),
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
			# title!("\n\nopt dispatch for JCC model, p = $(inputs.p), κ = $(inputs.κ), obj value = $(-round(objective_value(m)*1000, digits=3)) €", titlelocation=:center) #, fontsize=3)
			plot!(inputs.d, label = "forecast. demand", linewidth=2, color=:black)
			plot!(inputs.s_max, label = "forecast. total solar", linewidth=2, color=:orange)

		end 
	end 
	savefig("JCCModel_$(inputs.p *100)_$(inputs.κ).png")

	# Plot of reserves
	begin
		if termination_status(m) in [MOI.LOCALLY_SOLVED, MOI.OPTIMAL] 
			gr(size=(800,600))
			ticklabel = [1:length(T)+inputs.κ]
			batt = groupedbar([value.(m[:r_b_minus]) value.(m[:r_Δ])],
					bar_position = :stack,
					bar_width=0.7,
					xticks=(1:2:length(T)+inputs.κ),
					ylims=(0, 20),
					label=["battery discharge reserve" "diesel reserve"],
					legend=:topleft,
					top_margin = 10mm,
					xlabel = "\n\nTime in h",
					ylabel = "\n\nPower in kW",
					labelfontsize = 13,
					legendfontsize=12,
					tickfontsize = 11,
					color = [:magenta4 :thistle])
			# title!("\n\nopt reserves for JCC model, p = $(inputs.p), κ = $(inputs.κ), obj value = $(-round(objective_value(m)*1000, digits=3)) €", titlelocation=:center) #, fontsize=3)
		end 
	end
	savefig("Reserves_JCCModel_$(inputs.p *100)_$(inputs.κ).png")

	# Plot of battery profile
	begin
		if termination_status(m) in [MOI.LOCALLY_SOLVED, MOI.OPTIMAL] 
			gr(size=(900,700))
			ticklabel = [1:length(T)+inputs.κ]
			res = groupedbar([value.(m[:b_minus]) -value.(m[:b_plus])],
					bar_position = :stack,
					bar_width=0.7,
					xticks=(1:2:length(T)+inputs.κ),
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
			# title!("\n\nopt battery dispatch for JCC model, p = $(inputs.p), κ = $(inputs.κ), obj value = $(-round(objective_value(m)*1000, digits=3)) €", titlelocation=:center) #, fontsize=3)
			plot!(twinx(), value.(m[:SOC]), ylims=(15, 80), label = "SOC", linewidth=2, color=:black, ylabel="\n\nSOC in %\n\n", legendfontsize=12, tickfontsize = 14) #, legend=:bottomright)
		end 
	end
	savefig("Battery_JCCModel_$(inputs.p *100)_$(inputs.κ).png")

	begin
		Random.seed!(1234)
		μ_δ = vec([0 for i in 1:length(T)+outage])
		Σ_δ = cov(transpose(d_error[1:length(T)+outage,:]))
		σ_δ = diag(Σ_δ, 0).^0.5
		δ = []
		for i in 1:length(T)+outage
			append!(δ, rand(Normal(μ_δ[i], σ_δ[i]), 1))
		end

		μ_ζ = vec([0 for i in 1:length(T)+outage])
		Σ_ζ = cov(transpose(s_dry_error))
		σ_ζ = diag(Σ_ζ, 0).^0.5
		ζ = []
		for i in 1:length(T)+outage
			append!(ζ, rand(Normal(μ_ζ[i], σ_ζ[i]), 1))
		end

		γ = δ + inputs.d - inputs.s_max - ζ - value.(m[:Δ]) - value.(m[:g_plus]) - value.(m[:b_minus]) + value.(m[:b_plus]) + value.(m[:g_minus])
	end
	# Plot real demand
	begin
		if termination_status(m) in [MOI.LOCALLY_SOLVED, MOI.OPTIMAL] 
			gr(size=(800,400))
			ticklabel = [1:length(T)+inputs.κ]
			d = groupedbar([γ+ζ value.(m[:Δ])+inputs.s_max+value.(m[:b_minus])+value.(m[:g_plus])-value.(m[:g_minus])-value.(m[:b_plus])],
					bar_position = :stack,
					bar_width=0.7,
					xticks=(1:2:length(T)+inputs.κ),
					ylims=(-20,15),
					label=["instantan. dispatch (γ+ξ)" "planned dispatch"],
					legend=:bottomright,
					top_margin = 10mm,
					xlabel = "\nTime in h",
					ylabel = "\n\nPower in kW",
					labelfontsize = 13,
					legendfontsize=10,
					tickfontsize = 11,
					color= [:turquoise :lavender])
			# title!("\n\nopt dispatch for JCC model, p = $(inputs.p), κ = $(inputs.κ), obj value = $(-round(objective_value(m)*1000, digits=3)) €", titlelocation=:center) #, fontsize=3)
			plot!(inputs.d + δ, label = "real demand (d+δ)", linewidth=2, color=:black, bottom_margin = 20px)
			plot!(inputs.d, label = "forecasted demand (d)", linewidth=2, color=:darkorange)
			# plot!(inputs.s_max, label = "forecast. total solar", linewidth=2, color=:orange)
		end 
	end 
	savefig("JCCModel_realdemand_$(c_γ_plus)_$(inputs.p *100)_$(inputs.κ).png")
end