# ======================================================================
# Utils.jl
# Useful functions and structs
# ======================================================================

mutable struct param
	c_Δ
	c_b
	c_g_plus
	p_g_minus
	p_d
	η_plus
	η_minus
	κ
	Δ_max
	b_max
	# b_cap # TODO add this
	SOC_max
	SOC_min
	SOC_0
	p
	d
	s_max
	g_max
	ω
	Σ
	μ
	σ
	a
end

param(;c_Δ,
	c_b,
	c_g_plus,
	p_g_minus,
	p_d,
	η_plus,
	η_minus,
	κ,
	Δ_max,
	b_max,
	SOC_max,
	SOC_min,
	SOC_0,
	p,
	d,
	s_max,
	g_max,
	ω,
	Σ,
	μ,
	σ,
	a) = param(c_Δ,
		c_b,
		c_g_plus,
		p_g_minus,
		p_d,
		η_plus,
		η_minus,
		κ,
		Δ_max,
		b_max,
		SOC_max,
		SOC_min,
		SOC_0,
		p,
		d,
		s_max,
		g_max,
		ω,
		Σ,
		μ,
		σ,
		a
	)


	function dump_results(T, m, season, inputs)
	results_df = DataFrames.DataFrame(
			    Δ_t = Float64[],
				g_plus_t = Float64[],
				g_minus_t = Float64[],
				b_plus_t = Float64[],
				b_minus_t = Float64[],
				SOC_t = Float64[],
				s_t = Float64[],
				s_max_t = Float64[],
				d_t = Float64[],
				# ω_t = Float64[],
				objective = Float64[],
	)
	push!(results_df, (0, 0, 0, 0, 0, 0, 0, 0, 0, objective_value(m)))

	# TODO make this more efficient
	for i in 1:length(T)+inputs.κ
		push!(results_df, (value.(m[:Δ])[i], value.(m[:g_plus])[i], value.(m[:g_minus])[i], value.(m[:b_plus])[i], value.(m[:b_minus])[i], value.(m[:SOC])[i], value.(m[:s])[i], inputs.s_max[i], inputs.d[i], 0))
	end

	if isfile("Results_regular_season$(lpad(season, 1)).xlsx")
		rm("Results_regular_season$(lpad(season, 1)).xlsx")
	end
	XLSX.writetable("Results_regular_season$(lpad(season, 1)).xlsx", results = (collect(eachcol(results_df)), names(results_df)))
end

function dump_results_ind(T, m, season, inputs)
	results_df = DataFrames.DataFrame(
			    Δ_t = Float64[],
				g_plus_t = Float64[],
				g_minus_t = Float64[],
				b_plus_t = Float64[],
				b_minus_t = Float64[],
				r_Δ_t = Float64[],
				r_b_plus_t = Float64[],
				r_b_minus_t = Float64[],
				SOC_t = Float64[],
				s_t = Float64[],
				s_max_t = Float64[],
				d_t = Float64[],
				# ω_t = Float64[],
				objective = Float64[],
	)
	push!(results_df, (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, objective_value(m)))

	# TODO make this more efficient
	for i in 1:length(T)+inputs.κ
		push!(results_df, (value.(m[:Δ])[i], value.(m[:g_plus])[i], value.(m[:g_minus])[i], value.(m[:b_plus])[i], value.(m[:b_minus])[i], value.(m[:r_Δ])[i], value.(m[:r_b_plus])[i], value.(m[:r_b_minus])[i], value.(m[:SOC])[i], value.(m[:s])[i], inputs.s_max[i], inputs.d[i], 0))
	end

	if isfile("Results_individualCC_season$(lpad(season, 1))_prob$(lpad(inputs.p, 3)).xlsx")
		rm("Results_individualCC_season$(lpad(season, 1))_prob$(lpad(inputs.p, 3)).xlsx")
	end
	XLSX.writetable("Results_individualCC_season$(lpad(season, 1))_prob$(lpad(inputs.p, 3)).xlsx", results = (collect(eachcol(results_df)), names(results_df)))
end

function dump_results_joint(T, m, season, inputs)
	results_df = DataFrames.DataFrame(
			    Δ_t = Float64[],
				g_plus_t = Float64[],
				g_minus_t = Float64[],
				b_plus_t = Float64[],
				b_minus_t = Float64[],
				r_Δ_t = Float64[],
				r_b_plus_t = Float64[],
				r_b_minus_t = Float64[],
				SOC_t = Float64[],
				s_t = Float64[],
				s_max_t = Float64[],
				d_t = Float64[],
				# ω_t = Float64[],
				objective = Float64[],
	)
	push!(results_df, (0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, objective_value(m)))

	# TODO make this more efficient
	for i in 1:length(T)+inputs.κ
		push!(results_df, (value.(m[:Δ])[i], value.(m[:g_plus])[i], value.(m[:g_minus])[i], value.(m[:b_plus])[i], value.(m[:b_minus])[i], value.(m[:r_Δ])[i], value.(m[:r_b_plus])[i], value.(m[:r_b_minus])[i], value.(m[:SOC])[i], value.(m[:s])[i], inputs.s_max[i], inputs.d[i], 0))
	end
	if isfile("Results_JointCC_season$(lpad(season, 1))_kappa$(lpad(inputs.κ, 2, "0"))_prob$(lpad(inputs.p, 3)).xlsx")
		rm("Results_JointCC_season$(lpad(season, 1))_kappa$(lpad(inputs.κ, 2, "0"))_prob$(lpad(inputs.p, 3)).xlsx")
    end
	XLSX.writetable("Results_JointCC_season$(lpad(season, 1))_kappa$(lpad(inputs.κ, 2, "0"))_prob$(lpad(inputs.p, 3)).xlsx", results = (collect(eachcol(results_df)), names(results_df)))
end

function count_violations(m, a, f)
	counter = 0
	# r = value.(r_Δ) - value.(r_b_plus) + value.(r_b_minus)
	r = value.(r_Δ) + value.(r_b_minus)
	for i in 1:nsims
		if r[f] - a[f, i] < 0
			counter = counter +1
		end
	end
	return counter
end


function plot_results_individ(T, m, 𝛏, κ, nsims)
	x = value.(r_Δ) + value.(r_b_minus) - value.(g_plus) + value.(g_minus)
	fig_individ= plot()
	for i in 1:nsims
		plot!(fig_individ, x - 𝛏[:, i], )
	end
	plot!(fig_individ, xlabel = "\n\nTime in h\n\n", ylabel = "\n\nForecasting errors with reserves in kW", xticks = (1:length(T)+κ,[string(i) for i in 1:length(T)+κ]), ylims = (-15, 20), legend = false, framstyle = :zerolines, linewidth=2)

	fig_errors = plot()
	for i in 1:nsims
		plot!(fig_errors, - 𝛏[:, i],)
	end
	plot!(fig_errors, xlabel = "\n\nTime in h\n\n", ylabel = "\n\nForecasting errors without reserves in kW", xticks = (1:length(x)+κ,[string(i) for i in 1:length(x)+κ]), ylims = (-15, 20), legend = false, framstyle = :zerolines, linewidth=2)

	return fig_individ, fig_errors
end

function plot_results(T, m, 𝛏, f, κ, nsims)
	x = value.(m[:r_Δ]) + value.(m[:r_b_minus]) - value.(m[:g_plus]) + value.(m[:g_minus])

	fig_joint = plot()
	for i in 1:nsims
		plot!(fig_joint, x[f:f+κ] - 𝛏[f:f+κ, i], )
	end
	plot!(fig_joint, xlabel = "\nTime in h", ylabel = "\n\nForecasting errors with reserves in kW", xticks = (1:κ+1,[string(i) for i in f:f+κ]), ylims = (-15, 20), legend = false, framstyle = :zerolines, linewidth=2)

	fig_errors = plot()
	for i in 1:nsims
		plot!(fig_errors, - 𝛏[f:f+κ, i],)
	end
	plot!(fig_errors, xlabel = "\nTime in h", ylabel = "\n\nForecasting errors without reserves in kW", xticks = (1:κ+1, [string(i) for i in f:f+κ]), ylims = (-15, 20), legend = false, framstyle = :zerolines, linewidth=2)

	return fig_joint, fig_errors #, fig
end

# Create random sample
function create_sample(ξ, κ, nsims)
	Random.seed!(1234)
	γ = Array{Float64}(undef, nsims, 0)
	for i = 1:length(T)+κ 
		# γ = [γ rand(ξ[i] * σ[i] + μ[i], nsims)]
		γ = [γ rand(ξ[i], nsims)]
	end
	γ = transpose(γ)
	df = DataFrame(γ, :auto)
	if isfile("errors.xlsx")
		rm("errors.xlsx")
	end
	XLSX.writetable("errors.xlsx", collect(eachcol(df)), names(df))
	return γ
end

# Create random sample
function create_sample_joint(ξ, Σ, κ, nsims)
	Random.seed!(1234)
	γ = Array{Float64}(undef, nsims, 0)
	for i = 1:length(T)
		# γ = [γ rand(ξ[i] * σ[i] + μ[i], nsims)]
		# β = MvNormal([mean(ξ[x]) for x in 1:1+κ], Σ[1:1+κ, 1:1+κ])
		β = MvNormal([mean(ξ[x]) for x in i:i+κ], Σ[i:i+κ, i:i+κ])
		γ = [γ transpose(rand(β, nsims))]
	end
	γ = transpose(γ)
	df = DataFrame(γ, :auto)
	if isfile("errors_joint.xlsx")
		rm("errors_joint.xlsx")
	end
	XLSX.writetable("errors_joint.xlsx", collect(eachcol(df)), names(df))
	return γ
end

function Scenario(horizon::Vector, outage_duration::Int64, demand_scaling::Float64, season::Int64)
	return (horizon = horizon, outage_duration = outage_duration, demand_scaling = demand_scaling, season = season)
end
