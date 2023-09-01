# ======================================================================
# JointProb.jl
# This code is used to define a minigrid dispatch optimization model 
# with Joint Probabilistic (Chance) Constraints.
# ======================================================================

function define_model_JCC_Genz(scenario, inputs::param, m)
	T = scenario[1]
	outage = scenario[2]
	q = scenario[3]
	season = scenario[4]

	### Define positive decision variables
	@variable(m, 0. <= Δ[1:length(T)+inputs.κ] <= inputs.Δ_max)
	@variable(m, 0. <= b_plus[1:length(T)+inputs.κ] <= inputs.b_max)
	@variable(m, 0. <= b_minus[1:length(T)+inputs.κ] <= inputs.b_max)
	@variable(m, 0. <= g_plus[i = 1:length(T)+inputs.κ] <= inputs.g_max[i])
	@variable(m, 0. <= g_minus[i = 1:length(T)+inputs.κ] <= inputs.g_max[i])
	@variable(m, r_Δ[1:length(T)+inputs.κ] >= 0.)
	@variable(m, r_b_plus[1:length(T)+inputs.κ] >= 0.)
	@variable(m, r_b_minus[1:length(T)+inputs.κ] >= 0.)
	@variable(m, SOC[1:length(T)+inputs.κ] >= inputs.SOC_min)
	@variable(m, 0. <= s[i = 1:length(T)+inputs.κ] <= inputs.s_max[i])

	### add start values
	set_start_value.(Δ, map(+, [0 for i in 1:length(T)+outage], 
							repeat(round.([3.5278636760657256, 3.32526273300156, 2.682812817630818, 3.208789574358008, 3.1801556140373415, 3.0917609564919553, 3.133686815024245, 3.9728526512476923, 3.9600548440504044e-8, 3.9600536690735684e-8, 3.960060626876671e-8, 3.583714072067976e-8, 3.583705950958135e-8, 3.5836990127366804e-8, 3.5837061347719454e-8, 3.583703516732931e-8, 3.583713085137936e-8, 3.583701353788521e-8, 3.58370893219171e-8, 1.0873451813532386e-7, 1.1121211688676264e-7, 1.1121251807304423e-7, 1.1121161508454557e-7, 8.131060251430068e-8]), outer = 2)))	
	set_start_value.(b_plus, map(+, [0 for i in 1:length(T)+outage], 
							 repeat(round.([8.090909789445457e-8, 8.090909854343228e-8, 8.090910128868296e-8, 8.090909732987539e-8, 8.090909651418454e-8, 8.090909591023901e-8, 8.090909414759488e-8, 8.090907824350485e-8, 1.9512330708564727e-6, 9.675902311868854e-7, 7.16366316197616e-7, 4.44008240069916, 7.293494220461738, 6.404452980310228, 5.536464419961817, 4.843323271342555, 4.297056585905067, 3.4870353428502665, 2.2985034732678247, 5.4203361986463985e-8, 5.282722502394868e-8, 5.2827127930153155e-8, 5.2827364350427063e-8, 8.090907945350747e-8]), outer = 2)))	
	set_start_value.(b_minus, map(+, [0 for i in 1:length(T)+outage], repeat(round.([1.215206664904904, 1.0597854844951677, 0.7597712037118182, 1.1655806126384718, 1.3089434626544998, 1.4011751940955024, 1.724302628474902, 7.154709207769594, 8.532783902542512e-8, 9.022977312587494e-8, 9.391463597943235e-8, 8.090908649191675e-8, 8.090908119209515e-8, 8.090908462035346e-8, 8.090908642796437e-8, 8.090908787295802e-8, 8.09090874163301e-8, 8.090908720831687e-8, 8.090908753730865e-8, 2.6678128052224244, 6.953056207688357, 5.410364065162251, 7.247793305691659, 5.02524166389183]), outer = 2)))		
	set_start_value.(g_plus, map(+, [0 for i in 1:length(T)+outage], repeat(round.([4.647059072667199e-8, 5.48648673514613e-8, 6.619047697702985e-8, 8.230769654625505e-8, 8.230769618903475e-8, 8.230769523892873e-8, 8.230769472217081e-8, 8.230770720493615e-8, 12.907756721995352, 8.489212875147256, 1.1253130470602377, 5.900019340484925e-7, 5.900008439237207e-7, 5.900006094193698e-7, 5.900005076542211e-7, 5.900004588850108e-7, 5.90000434317633e-7, 5.900004528132629e-7, 1.8571428999770516e-8, 4.148139993024351e-8, 4.240174788930317e-8, 4.240181615937557e-8, 4.240165172630243e-8, 3.166667027081109e-8]), outer = 2)))		
	set_start_value.(g_minus, map(+, [0 for i in 1:length(T)+outage], repeat(round.([2.6585364801053256e-8, 2.614457753627444e-8, 2.5714285538381436e-8, 2.5294117029576914e-8, 2.5294117081662283e-8, 2.529411722011098e-8, 2.529411729528199e-8, 2.5294115477129813e-8, 5.900006673831758e-7, 5.900007452189407e-7, 5.900010062592879e-7, 0.36406341118193336, 3.107698774425298, 8.110254464163182, 11.266757583968941, 12.040556774264815, 11.31496342707562, 6.724678339949843, 3.032617789635995, 5.420336028040823e-8, 5.2827222695055164e-8, 3.9792468415854604e-8, 3.979261678080308e-8, 5.5934056879441005e-8]), outer = 2)))	
	set_start_value.(r_Δ, map(+, [0 for i in 1:length(T)+outage], repeat(round.([8.035614513546174e-7, 3.9678032433936406e-7, 2.61187110969062e-7, 1.9338925818922668e-7, 1.9338961244409089e-7, 1.9338990575767337e-7, 1.9338986523012482e-7, 1.9338801177291754e-7, 3.77943345886614e-7, 3.7794406239321077e-7, 3.7793940644610454e-7, 3.796426894920193e-7, 3.796506119948257e-7, 3.7965401792323563e-7, 3.7965098080625117e-7, 3.796478913510762e-7, 3.7964305541334895e-7, 3.7965165733461276e-7, 3.796439556075238e-7, 5.798248857141978e-6, 3.6366531733476952, 0.46848294019493225, 1.6402744296590002, 3.60148404522889e-7]), outer = 2)))			
	set_start_value.(r_b_plus, map(+, [0 for i in 1:length(T)+outage], repeat(round.([10.000 for i in 1:length(T)]), outer = 2)))	
	set_start_value.(r_b_minus, map(+, [0 for i in 1:length(T)+outage], repeat(round.([2.1655917311523143, 2.87626444868623, 3.958760722483419, 3.4831456824139764, 3.2382491151106194, 3.0292595301621676, 2.6676148774404895, 2.2585023854757607, 3.056611854410409, 4.357863303760356, 5.398135073501355, 5.171421004969157, 5.606716575123699, 6.105230050986148, 6.856825784576239, 6.872147784941361, 7.834833156039777, 8.358622439825707, 8.064368036819623, 7.332187076222642, 3.0469436889218096, 4.589635831448936, 2.752206590917041, 4.268721587341084]), outer = 2)))	
	set_start_value.(SOC, map(+, [0 for i in 1:length(T)+outage], repeat(round.([33.84555374520398, 32.83875761179722, 32.11697504513464, 31.00967353999173, 29.766177327333597, 28.435060969806507, 26.796973549618993, 19.9999998791015, 20.00000165171147, 20.000002485203904, 20.000003076533, 24.218081280333575, 31.146900712908597, 37.23113096733968, 42.49077208943977, 47.091929120351566, 51.17413280009774, 54.48681629894186, 56.67039452168266, 54.13597240821455, 47.53056906109647, 42.3907232493781, 35.505319659157024, 30.73134015532341]), outer = 2)))	
	set_start_value.(s, inputs.s_max)

	### Define the objective function
	@objective(
		m, 
		Min, 
		(sum(inputs.c_Δ .* Δ
			+ inputs.c_b * (b_plus + b_minus)
			- inputs.p_d * inputs.d)
			+ (sum([sum([inputs.c_g_plus[t] * g_plus[t] - inputs.p_g_minus[t] * g_minus[t] 
									   for t in 1:length(T)+inputs.κ if t ∉ τ:τ+inputs.κ])
									   for τ in 1:length(T)]) 
			    + (sum([sum([inputs.c_Δ .* r_Δ[t] + inputs.c_b .* r_b_minus[t] for t in τ:τ+inputs.κ]) for τ in eachindex(T)]))) * (inputs.ω  / (length(T)+inputs.κ)) 
			+ sum([inputs.c_g_plus[t] * g_plus[t] - inputs.p_g_minus[t] * g_minus[t] for t in 1:length(T)+inputs.κ]) * (1-inputs.ω)
		) * 0.001 # scaling objective function to help solver
		)

	### Define the constraints
	@constraint(m, Δ + r_Δ .- inputs.Δ_max .<= 0.00)
	@constraint(m, b_minus + r_b_minus .- inputs.b_max .<= 0.00)
	@constraint(m, [t = 2:length(T)+inputs.κ], 
					SOC[t] == SOC[t-1] + (inputs.η_plus * b_plus[t] - inputs.η_minus * b_minus[t]) * Δt / b_cap
				)
	@constraint(m, 
				SOC[1] == inputs.SOC_0 + (inputs.η_plus * b_plus[1] - inputs.η_minus * b_minus[1]) * Δt / b_cap
				)
	@constraint(m, [t = 2:length(T)+inputs.κ], 
					inputs.SOC_min <= SOC[t-1] + (inputs.η_plus * (b_plus[t]) + inputs.η_minus * (b_minus[t] + 
					r_b_minus[t])) * Δt / b_cap
				)
	@constraint(m, [t = 2:length(T)+inputs.κ], 
				SOC[t-1] + (inputs.η_plus * (b_plus[t]) + inputs.η_minus * (b_minus[t] + 
				r_b_minus[t])) * Δt / b_cap .- inputs.SOC_max .<= 0.00
				)
	@constraint(m, inputs.SOC_min <= inputs.SOC_0 + (inputs.η_plus * (b_plus[1]) +
				   inputs.η_minus * (b_minus[1] + r_b_minus[1])) * Δt / b_cap
				)
	@constraint(m, inputs.SOC_0 + (inputs.η_plus * (b_plus[1]) +
				inputs.η_minus * (b_minus[1] + r_b_minus[1])) * Δt / b_cap - inputs.SOC_max <= 0.00
			 	)
	@constraint(m, s + Δ + g_plus + b_minus - b_plus - g_minus - inputs.d .== 0)
	@constraint(m, SOC[length(T)] == inputs.SOC_0)

	function distribution(j, inputs)
		norm_pdf(k) = exp(-(k^2) / 2) / (2 * pi)^0.5
		f(x...) = mvnormcdf(inputs.Σ[j:j+inputs.κ, j:j+inputs.κ], inputs.a[j:j+inputs.κ]-inputs.μ[j:j+inputs.κ], vec([x[i] for i in j:j+inputs.κ])-inputs.μ[j:j+inputs.κ], m = 5000, rng = MersenneTwister(1234))[1]
		function ∇f(g::AbstractVector{T}, x::T...) where {T}
			μ_j = inputs.μ[j:j+inputs.κ]
			σ_j = inputs.σ[j:j+inputs.κ]
			Σ_j = inputs.Σ[j:j+inputs.κ, j:j+inputs.κ]
			a_j = inputs.a[j:j+inputs.κ]
			for i in eachindex(x)
                if i in j:j+inputs.κ
                    Σ_new = Σ_j - inv(Σ_j[i-j+1, i-j+1]) * Σ_j[i-j+1, :] * transpose(Σ_j[i-j+1,:])
                    μ_new = μ_j + inv(Σ_j[i-j+1, i-j+1]) * (x[i] - μ_j[i-j+1]) * Σ_j[i-j+1, :]
                    g[i] = norm_pdf((x[i] - μ_j[i-j+1])/σ_j[i-j+1]) / σ_j[i-j+1] * mvnormcdf(Σ_new[1:end .!= i-j+1, 
																								   1:end .!= i-j+1], a_j[1:end .!=i-j+1]-μ_new[1:end .!=i-j+1], 
																								   vec([x[k] for k in j:j+inputs.κ if k != i])-μ_new[1:end .!=i-j+1], 
																								   m = 5000, rng = MersenneTwister(1234))[1]  
                else 
                    g[i] = 0.00000
                end
            end
            return
		end
		return f, ∇f
	end

	x = r_b_minus + r_Δ - g_plus + g_minus

	for j in eachindex(T)
		register(m, Symbol("mvncdf_$j"), length(x), distribution(j, inputs)[1], distribution(j, inputs)[2])
		add_nonlinear_constraint(m, :($(Symbol("mvncdf_$j"))($(x...)) >= $(inputs.p)))
	end

	return m
end


function define_model_JCC_SRD(scenario, inputs::param, m)

	begin
		dir = @__DIR__
		include("$dir/SRD/chi_gamma_cdf.jl")
		include("$dir/SRD/chi_gamma_pdf.jl")
		include("$dir/SRD/SRDFunction.jl")
		include("$dir/SRD/SetProblemData.jl")
	end

	T = scenario[1]
	outage = scenario[2]
	q = scenario[3]
	season = scenario[4]

	#### Define positive decision variables
	@variable(m, 0. <= Δ[1:length(T)+inputs.κ] <= inputs.Δ_max)
	@variable(m, 0. <= b_plus[1:length(T)+inputs.κ] <= inputs.b_max)
	@variable(m, 0. <= b_minus[1:length(T)+inputs.κ] <= inputs.b_max)
	@variable(m, 0. <= g_plus[i = 1:length(T)+inputs.κ] <= inputs.g_max[i])
	@variable(m, 0. <= g_minus[i = 1:length(T)+inputs.κ] <= inputs.g_max[i])
	@variable(m, r_Δ[1:length(T)+inputs.κ] >= 0.)
	@variable(m, r_b_plus[1:length(T)+inputs.κ] >= 0.)
	@variable(m, r_b_minus[1:length(T)+inputs.κ] >= 0.)
	@variable(m, SOC[1:length(T)+inputs.κ] >= inputs.SOC_min)
	@variable(m, 0. <= s[i = 1:length(T)+inputs.κ] <= inputs.s_max[i])

	### add start values
	set_start_value.(Δ, map(+, [0 for i in 1:length(T)+outage], 
							repeat(round.([3.5278636760657256, 3.32526273300156, 2.682812817630818, 3.208789574358008, 3.1801556140373415, 3.0917609564919553, 3.133686815024245, 3.9728526512476923, 3.9600548440504044e-8, 3.9600536690735684e-8, 3.960060626876671e-8, 3.583714072067976e-8, 3.583705950958135e-8, 3.5836990127366804e-8, 3.5837061347719454e-8, 3.583703516732931e-8, 3.583713085137936e-8, 3.583701353788521e-8, 3.58370893219171e-8, 1.0873451813532386e-7, 1.1121211688676264e-7, 1.1121251807304423e-7, 1.1121161508454557e-7, 8.131060251430068e-8]), outer = 2)))	
	set_start_value.(b_plus, map(+, [0 for i in 1:length(T)+outage], 
							 repeat(round.([8.090909789445457e-8, 8.090909854343228e-8, 8.090910128868296e-8, 8.090909732987539e-8, 8.090909651418454e-8, 8.090909591023901e-8, 8.090909414759488e-8, 8.090907824350485e-8, 1.9512330708564727e-6, 9.675902311868854e-7, 7.16366316197616e-7, 4.44008240069916, 7.293494220461738, 6.404452980310228, 5.536464419961817, 4.843323271342555, 4.297056585905067, 3.4870353428502665, 2.2985034732678247, 5.4203361986463985e-8, 5.282722502394868e-8, 5.2827127930153155e-8, 5.2827364350427063e-8, 8.090907945350747e-8]), outer = 2)))	
	set_start_value.(b_minus, map(+, [0 for i in 1:length(T)+outage], repeat(round.([1.215206664904904, 1.0597854844951677, 0.7597712037118182, 1.1655806126384718, 1.3089434626544998, 1.4011751940955024, 1.724302628474902, 7.154709207769594, 8.532783902542512e-8, 9.022977312587494e-8, 9.391463597943235e-8, 8.090908649191675e-8, 8.090908119209515e-8, 8.090908462035346e-8, 8.090908642796437e-8, 8.090908787295802e-8, 8.09090874163301e-8, 8.090908720831687e-8, 8.090908753730865e-8, 2.6678128052224244, 6.953056207688357, 5.410364065162251, 7.247793305691659, 5.02524166389183]), outer = 2)))		
	set_start_value.(g_plus, map(+, [0 for i in 1:length(T)+outage], repeat(round.([4.647059072667199e-8, 5.48648673514613e-8, 6.619047697702985e-8, 8.230769654625505e-8, 8.230769618903475e-8, 8.230769523892873e-8, 8.230769472217081e-8, 8.230770720493615e-8, 12.907756721995352, 8.489212875147256, 1.1253130470602377, 5.900019340484925e-7, 5.900008439237207e-7, 5.900006094193698e-7, 5.900005076542211e-7, 5.900004588850108e-7, 5.90000434317633e-7, 5.900004528132629e-7, 1.8571428999770516e-8, 4.148139993024351e-8, 4.240174788930317e-8, 4.240181615937557e-8, 4.240165172630243e-8, 3.166667027081109e-8]), outer = 2)))		
	set_start_value.(g_minus, map(+, [0 for i in 1:length(T)+outage], repeat(round.([2.6585364801053256e-8, 2.614457753627444e-8, 2.5714285538381436e-8, 2.5294117029576914e-8, 2.5294117081662283e-8, 2.529411722011098e-8, 2.529411729528199e-8, 2.5294115477129813e-8, 5.900006673831758e-7, 5.900007452189407e-7, 5.900010062592879e-7, 0.36406341118193336, 3.107698774425298, 8.110254464163182, 11.266757583968941, 12.040556774264815, 11.31496342707562, 6.724678339949843, 3.032617789635995, 5.420336028040823e-8, 5.2827222695055164e-8, 3.9792468415854604e-8, 3.979261678080308e-8, 5.5934056879441005e-8]), outer = 2)))	
	set_start_value.(r_Δ, map(+, [0 for i in 1:length(T)+outage], repeat(round.([8.035614513546174e-7, 3.9678032433936406e-7, 2.61187110969062e-7, 1.9338925818922668e-7, 1.9338961244409089e-7, 1.9338990575767337e-7, 1.9338986523012482e-7, 1.9338801177291754e-7, 3.77943345886614e-7, 3.7794406239321077e-7, 3.7793940644610454e-7, 3.796426894920193e-7, 3.796506119948257e-7, 3.7965401792323563e-7, 3.7965098080625117e-7, 3.796478913510762e-7, 3.7964305541334895e-7, 3.7965165733461276e-7, 3.796439556075238e-7, 5.798248857141978e-6, 3.6366531733476952, 0.46848294019493225, 1.6402744296590002, 3.60148404522889e-7]), outer = 2)))			
	set_start_value.(r_b_plus, map(+, [0 for i in 1:length(T)+outage], repeat(round.([10.000 for i in 1:length(T)]), outer = 2)))	
	set_start_value.(r_b_minus, map(+, [0 for i in 1:length(T)+outage], repeat(round.([2.1655917311523143, 2.87626444868623, 3.958760722483419, 3.4831456824139764, 3.2382491151106194, 3.0292595301621676, 2.6676148774404895, 2.2585023854757607, 3.056611854410409, 4.357863303760356, 5.398135073501355, 5.171421004969157, 5.606716575123699, 6.105230050986148, 6.856825784576239, 6.872147784941361, 7.834833156039777, 8.358622439825707, 8.064368036819623, 7.332187076222642, 3.0469436889218096, 4.589635831448936, 2.752206590917041, 4.268721587341084]), outer = 2)))	
	set_start_value.(SOC, map(+, [0 for i in 1:length(T)+outage], repeat(round.([33.84555374520398, 32.83875761179722, 32.11697504513464, 31.00967353999173, 29.766177327333597, 28.435060969806507, 26.796973549618993, 19.9999998791015, 20.00000165171147, 20.000002485203904, 20.000003076533, 24.218081280333575, 31.146900712908597, 37.23113096733968, 42.49077208943977, 47.091929120351566, 51.17413280009774, 54.48681629894186, 56.67039452168266, 54.13597240821455, 47.53056906109647, 42.3907232493781, 35.505319659157024, 30.73134015532341]), outer = 2)))	
	set_start_value.(s, inputs.s_max)

	### Define the objective function
	@objective(
		m, 
		Min, 
		(sum(inputs.c_Δ .* Δ
			+ inputs.c_b * (b_plus + b_minus)
			- inputs.p_d * inputs.d)
			+ (length(T) == 1 ? 0 : sum([sum([inputs.c_g_plus[t] * g_plus[t] - inputs.p_g_minus[t] * g_minus[t] 
									   for t in 1:length(T)+inputs.κ if t ∉ τ:τ+inputs.κ])
									   for τ in 1:length(T)]) * (inputs.ω  / (length(T)+inputs.κ)))
			+ (sum([sum([inputs.c_Δ .* r_Δ[t] + inputs.c_b .* r_b_minus[t] for t in τ:τ+inputs.κ]) for τ in eachindex(T)])) * (inputs.ω  / (length(T)+inputs.κ)) 
		) * 0.001 # scaling objective function to help solver
		)

	### Define the constraints
	@constraint(m, Δ + r_Δ .- inputs.Δ_max .<= 0.00)
	@constraint(m, b_minus + r_b_minus .- inputs.b_max .<= 0.00)
	@constraint(m, [t = 2:length(T)+inputs.κ], 
					SOC[t] == SOC[t-1] + (inputs.η_plus * b_plus[t] - inputs.η_minus * b_minus[t]) * Δt / b_cap
				)
	@constraint(m, 
				SOC[1] == inputs.SOC_0 + (inputs.η_plus * b_plus[1] - inputs.η_minus * b_minus[1]) * Δt / b_cap
				)
	@constraint(m, [t = 2:length(T)+inputs.κ], 
					inputs.SOC_min <= SOC[t-1] + (inputs.η_plus * (b_plus[t]) + inputs.η_minus * (b_minus[t] + 
					r_b_minus[t])) * Δt / b_cap
				)
	@constraint(m, [t = 2:length(T)+inputs.κ], 
				SOC[t-1] + (inputs.η_plus * (b_plus[t]) + inputs.η_minus * (b_minus[t] + 
				r_b_minus[t])) * Δt / b_cap .- inputs.SOC_max .<= 0.00
				)
	@constraint(m, inputs.SOC_min <= inputs.SOC_0 + (inputs.η_plus * (b_plus[1]) +
				   inputs.η_minus * (b_minus[1] + r_b_minus[1])) * Δt / b_cap
				)
	@constraint(m, inputs.SOC_0 + (inputs.η_plus * (b_plus[1]) +
				inputs.η_minus * (b_minus[1] + r_b_minus[1])) * Δt / b_cap - inputs.SOC_max <= 0.00
			 	)
	@constraint(m, s + Δ + g_plus + b_minus - b_plus - g_minus - inputs.d .== 0)

	function create_cons(w, inputs)
		k = (1+inputs.κ) * (w-1) + 1
		f(x...) = ProbFunction(vec([x[r] for r in w:w+inputs.κ]), SampleOnSphere[:, k:k+inputs.κ], inputs.μ[w:w+inputs.κ])[1]
		function ∇f(g::AbstractVector{T}, x::T...) where {T}
			for i in eachindex(x)
				if i in w:w+inputs.κ
					g[i] = ProbFunction(vec([x[r] for r in w:w+inputs.κ]), SampleOnSphere[:, k:k+inputs.κ], inputs.μ[w:w+inputs.κ])[2][i-w+1]
				else 
					g[i] = 0.000
				end 
			end
			return
		end
		return f, ∇f
	end 
	
	x = r_b_minus + r_Δ - g_plus + g_minus

	for j in eachindex(T)
		register(m, Symbol("srd_prob_$j"), length(x), create_cons(j, inputs)[1], create_cons(j, inputs)[2])
		add_nonlinear_constraint(m, :($(Symbol("srd_prob_$j"))($(x...)) >= $(inputs.p)))
	end

	return m
end

function define_p_max_problem(scenario, inputs::param, m)
	T = scenario[1]
	outage = scenario[2]
	q = scenario[3]
	season = scenario[4]

	### Define positive decision variables
	@variable(m, 0. <= Δ[1:length(T)+inputs.κ] <= inputs.Δ_max)
	@variable(m, 0. <= b_plus[1:length(T)+inputs.κ] <= inputs.b_max)
	@variable(m, 0. <= b_minus[1:length(T)+inputs.κ] <= inputs.b_max)
	@variable(m, 0. <= g_plus[i = 1:length(T)+inputs.κ] <= inputs.g_max[i])
	@variable(m, 0. <= g_minus[i = 1:length(T)+inputs.κ] <= inputs.g_max[i])
	@variable(m, r_Δ[1:length(T)+inputs.κ] >= 0.)
	@variable(m, r_b_plus[1:length(T)+inputs.κ] >= 0.)
	@variable(m, r_b_minus[1:length(T)+inputs.κ] >= 0.)
	@variable(m, SOC[1:length(T)+inputs.κ] >= inputs.SOC_min)
	@variable(m, 0. <= s[i = 1:length(T)+inputs.κ] <= inputs.s_max[i])
	@variable(m, 0. <= p <= 1., start = 0.9)

	### add start values
	set_start_value.(Δ, map(+, [0 for i in 1:length(T)+outage], 
							repeat(round.([3.5278636760657256, 3.32526273300156, 2.682812817630818, 3.208789574358008, 3.1801556140373415, 3.0917609564919553, 3.133686815024245, 3.9728526512476923, 3.9600548440504044e-8, 3.9600536690735684e-8, 3.960060626876671e-8, 3.583714072067976e-8, 3.583705950958135e-8, 3.5836990127366804e-8, 3.5837061347719454e-8, 3.583703516732931e-8, 3.583713085137936e-8, 3.583701353788521e-8, 3.58370893219171e-8, 1.0873451813532386e-7, 1.1121211688676264e-7, 1.1121251807304423e-7, 1.1121161508454557e-7, 8.131060251430068e-8]), outer = 2)))	
	set_start_value.(b_plus, map(+, [0 for i in 1:length(T)+outage], 
							 repeat(round.([8.090909789445457e-8, 8.090909854343228e-8, 8.090910128868296e-8, 8.090909732987539e-8, 8.090909651418454e-8, 8.090909591023901e-8, 8.090909414759488e-8, 8.090907824350485e-8, 1.9512330708564727e-6, 9.675902311868854e-7, 7.16366316197616e-7, 4.44008240069916, 7.293494220461738, 6.404452980310228, 5.536464419961817, 4.843323271342555, 4.297056585905067, 3.4870353428502665, 2.2985034732678247, 5.4203361986463985e-8, 5.282722502394868e-8, 5.2827127930153155e-8, 5.2827364350427063e-8, 8.090907945350747e-8]), outer = 2)))	
	set_start_value.(b_minus, map(+, [0 for i in 1:length(T)+outage], repeat(round.([1.215206664904904, 1.0597854844951677, 0.7597712037118182, 1.1655806126384718, 1.3089434626544998, 1.4011751940955024, 1.724302628474902, 7.154709207769594, 8.532783902542512e-8, 9.022977312587494e-8, 9.391463597943235e-8, 8.090908649191675e-8, 8.090908119209515e-8, 8.090908462035346e-8, 8.090908642796437e-8, 8.090908787295802e-8, 8.09090874163301e-8, 8.090908720831687e-8, 8.090908753730865e-8, 2.6678128052224244, 6.953056207688357, 5.410364065162251, 7.247793305691659, 5.02524166389183]), outer = 2)))		
	set_start_value.(g_plus, map(+, [0 for i in 1:length(T)+outage], repeat(round.([4.647059072667199e-8, 5.48648673514613e-8, 6.619047697702985e-8, 8.230769654625505e-8, 8.230769618903475e-8, 8.230769523892873e-8, 8.230769472217081e-8, 8.230770720493615e-8, 12.907756721995352, 8.489212875147256, 1.1253130470602377, 5.900019340484925e-7, 5.900008439237207e-7, 5.900006094193698e-7, 5.900005076542211e-7, 5.900004588850108e-7, 5.90000434317633e-7, 5.900004528132629e-7, 1.8571428999770516e-8, 4.148139993024351e-8, 4.240174788930317e-8, 4.240181615937557e-8, 4.240165172630243e-8, 3.166667027081109e-8]), outer = 2)))		
	set_start_value.(g_minus, map(+, [0 for i in 1:length(T)+outage], repeat(round.([2.6585364801053256e-8, 2.614457753627444e-8, 2.5714285538381436e-8, 2.5294117029576914e-8, 2.5294117081662283e-8, 2.529411722011098e-8, 2.529411729528199e-8, 2.5294115477129813e-8, 5.900006673831758e-7, 5.900007452189407e-7, 5.900010062592879e-7, 0.36406341118193336, 3.107698774425298, 8.110254464163182, 11.266757583968941, 12.040556774264815, 11.31496342707562, 6.724678339949843, 3.032617789635995, 5.420336028040823e-8, 5.2827222695055164e-8, 3.9792468415854604e-8, 3.979261678080308e-8, 5.5934056879441005e-8]), outer = 2)))	
	set_start_value.(r_Δ, map(+, [0 for i in 1:length(T)+outage], repeat(round.([8.035614513546174e-7, 3.9678032433936406e-7, 2.61187110969062e-7, 1.9338925818922668e-7, 1.9338961244409089e-7, 1.9338990575767337e-7, 1.9338986523012482e-7, 1.9338801177291754e-7, 3.77943345886614e-7, 3.7794406239321077e-7, 3.7793940644610454e-7, 3.796426894920193e-7, 3.796506119948257e-7, 3.7965401792323563e-7, 3.7965098080625117e-7, 3.796478913510762e-7, 3.7964305541334895e-7, 3.7965165733461276e-7, 3.796439556075238e-7, 5.798248857141978e-6, 3.6366531733476952, 0.46848294019493225, 1.6402744296590002, 3.60148404522889e-7]), outer = 2)))			
	set_start_value.(r_b_plus, map(+, [0 for i in 1:length(T)+outage], repeat(round.([10.000 for i in 1:length(T)]), outer = 2)))	
	set_start_value.(r_b_minus, map(+, [0 for i in 1:length(T)+outage], repeat(round.([2.1655917311523143, 2.87626444868623, 3.958760722483419, 3.4831456824139764, 3.2382491151106194, 3.0292595301621676, 2.6676148774404895, 2.2585023854757607, 3.056611854410409, 4.357863303760356, 5.398135073501355, 5.171421004969157, 5.606716575123699, 6.105230050986148, 6.856825784576239, 6.872147784941361, 7.834833156039777, 8.358622439825707, 8.064368036819623, 7.332187076222642, 3.0469436889218096, 4.589635831448936, 2.752206590917041, 4.268721587341084]), outer = 2)))	
	set_start_value.(SOC, map(+, [0 for i in 1:length(T)+outage], repeat(round.([33.84555374520398, 32.83875761179722, 32.11697504513464, 31.00967353999173, 29.766177327333597, 28.435060969806507, 26.796973549618993, 19.9999998791015, 20.00000165171147, 20.000002485203904, 20.000003076533, 24.218081280333575, 31.146900712908597, 37.23113096733968, 42.49077208943977, 47.091929120351566, 51.17413280009774, 54.48681629894186, 56.67039452168266, 54.13597240821455, 47.53056906109647, 42.3907232493781, 35.505319659157024, 30.73134015532341]), outer = 2)))	
	set_start_value.(s, inputs.s_max)


	### Define the objective function
	@objective(
		m, 
		Max, 
		p*100 # scaling objective function to help solver
		)

	### Define the constraints
	@constraint(m, Δ + r_Δ .- inputs.Δ_max .<= 0.00)
	@constraint(m, b_minus + r_b_minus .- inputs.b_max .<= 0.00)
	@constraint(m, [t = 2:length(T)+inputs.κ], 
					SOC[t] == SOC[t-1] + (inputs.η_plus * b_plus[t] - inputs.η_minus * b_minus[t]) * Δt / b_cap
				)
	@constraint(m, 
				SOC[1] == inputs.SOC_0 + (inputs.η_plus * b_plus[1] - inputs.η_minus * b_minus[1]) * Δt / b_cap
				)
	@constraint(m, [t = 2:length(T)+inputs.κ], 
					inputs.SOC_min <= SOC[t-1] + (inputs.η_plus * (b_plus[t]) + inputs.η_minus * (b_minus[t] + 
					r_b_minus[t])) * Δt / b_cap
				)
	@constraint(m, [t = 2:length(T)+inputs.κ], 
				SOC[t-1] + (inputs.η_plus * (b_plus[t]) + inputs.η_minus * (b_minus[t] + 
				r_b_minus[t])) * Δt / b_cap .- inputs.SOC_max .<= 0.00
				)
	@constraint(m, inputs.SOC_min <= inputs.SOC_0 + (inputs.η_plus * (b_plus[1]) +
				   inputs.η_minus * (b_minus[1] + r_b_minus[1])) * Δt / b_cap
				)
	@constraint(m, inputs.SOC_0 + (inputs.η_plus * (b_plus[1]) +
				inputs.η_minus * (b_minus[1] + r_b_minus[1])) * Δt / b_cap - inputs.SOC_max <= 0.00
			 	)
	@constraint(m, s + Δ + g_plus + b_minus - b_plus - g_minus - inputs.d .== 0)
	@constraint(m, SOC[length(T)] == inputs.SOC_0)

	function distribution(j, inputs)
		norm_pdf(k) = exp(-(k^2) / 2) / (2 * pi)^0.5
		f(x...) = mvnormcdf(inputs.Σ[j:j+inputs.κ, j:j+inputs.κ], inputs.a[j:j+inputs.κ]-inputs.μ[j:j+inputs.κ], vec([x[i] for i in j:j+inputs.κ])-inputs.μ[j:j+inputs.κ], m = 5000, rng = MersenneTwister(1234))[1]
		function ∇f(g::AbstractVector{T}, x::T...) where {T}
			μ_j = inputs.μ[j:j+inputs.κ]
			σ_j = inputs.σ[j:j+inputs.κ]
			Σ_j = inputs.Σ[j:j+inputs.κ, j:j+inputs.κ]
			a_j = inputs.a[j:j+inputs.κ]
			for i in eachindex(x)
                if i in j:j+inputs.κ
                    Σ_new = Σ_j - inv(Σ_j[i-j+1, i-j+1]) * Σ_j[i-j+1, :] * transpose(Σ_j[i-j+1,:])
                    μ_new = μ_j + inv(Σ_j[i-j+1, i-j+1]) * (x[i] - μ_j[i-j+1]) * Σ_j[i-j+1, :]
                    g[i] = norm_pdf((x[i] - μ_j[i-j+1])/σ_j[i-j+1]) / σ_j[i-j+1] * mvnormcdf(Σ_new[1:end .!= i-j+1, 
																								   1:end .!= i-j+1], a_j[1:end .!=i-j+1]-μ_new[1:end .!=i-j+1], 
																								   vec([x[k] for k in j:j+inputs.κ if k != i])-μ_new[1:end .!=i-j+1], 
																								   m = 5000, rng = MersenneTwister(1234))[1]  
                else 
                    g[i] = 0.00000
                end
            end
            return
		end
		return f, ∇f
	end

	x = r_b_minus + r_Δ - g_plus + g_minus

	for j in eachindex(T)
		register(m, Symbol("mvncdf_$j"), length(x), distribution(j, inputs)[1], distribution(j, inputs)[2])
		add_nonlinear_constraint(m, :($(Symbol("mvncdf_$j"))($(x...)) >= $(p)))
	end

	return m
end