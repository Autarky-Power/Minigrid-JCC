## Define input parameters

include("$dir/InputDistributions.jl")

# Define parameters for JCC case
inputs = param(
	c_Δ = 0.35,
	c_b = 0.0055,
	c_g_plus = map(+, [0 for i in 1:length(T) + outage], 
			   repeat([0.55,0.55,0.55,0.55, 0.55, 0.55,
			   0.55, 0.55, 0.15, 0.15, 0.15, 0.15,
			   0.15, 0.15, 0.15, 0.15, 0.15, 0.15,
			   0.55, 0.55, 0.55, 0.55, 0.55, 0.55], outer = 2)),
	p_g_minus = map(+, [0 for i in 1:length(T) + outage], 
				repeat([0.08, 0.08, 0.08, 0.08, 0.08, 0.08,
				0.08, 0.08, 0.13, 0.13, 0.13, 0.13,
				0.13, 0.13, 0.13, 0.13, 0.13, 0.13,
				0.13, 0.13, 0.13, 0.08, 0.08, 0.08], outer = 2)),
	p_d = 0.55,
	η_plus = 0.95,
	η_minus = 0.95,
	κ = outage,
	Δ_max = 5,
	b_max = 10,  # TODO check this
	SOC_max = 90, # TODO 90
	SOC_min = 20, # TODO 20
	SOC_0 = 35, # TODO 35
	d = q * Float64.(d[:,end]),
	s_max = Float64.(s_dry[:,end]), # solar dry season
	g_max = [100 for i in 1:length(T)+outage],
	p = prob,
	ω = 0.9,
	Σ = Σ,
	μ = μ,
	σ = σ,
	a = a
	)

## Define parameters for ICC case
inputs_ind = param(
	c_Δ = 0.35,
	c_b = 0.0055,
	c_g_plus = map(+, [0 for i in 1:length(T) + outage], 
			   repeat([0.55,0.55,0.55,0.55, 0.55, 0.55,
			   0.55, 0.55, 0.15, 0.15, 0.15, 0.15,
			   0.15, 0.15, 0.15, 0.15, 0.15, 0.15,
			   0.55, 0.55, 0.55, 0.55, 0.55, 0.55], outer = 2)),
	p_g_minus = map(+, [0 for i in 1:length(T) + outage], 
				repeat([0.08, 0.08, 0.08, 0.08, 0.08, 0.08,
				0.08, 0.08, 0.13, 0.13, 0.13, 0.13,
				0.13, 0.13, 0.13, 0.13, 0.13, 0.13,
				0.13, 0.13, 0.13, 0.08, 0.08, 0.08], outer = 2)),
	p_d = 0.55,
	η_plus = 0.95,
	η_minus = 0.95,
	κ = outage,
	Δ_max = 5, 
	b_max = 10,
	SOC_max = 90,
	SOC_min = 20,
	SOC_0 = 35,
	d = q * Float64.(d[:, end]),
	s_max = Float64.(s_dry[:,end]), # solar dry season
	g_max = [100 for i in 1:length(T)+outage],
	p = prob,
	ω = .9,
	Σ = Σ,
	μ = μ,
	σ = σ,
	a = a
	)

c_γ_plus = 0.3
Δt = 1
# Battery capacity 100 kWh, divided by 100 to account for percentage values in kWh:
b_cap = 100/100
