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
	@variable(m, E[1:length(T)+inputs.κ] >= 0., start = 0.00)

	### add start values
	set_start_value.(Δ, map(+, [0 for i in 1:length(T)+outage], 
		repeat(round.([4.999999875405635, 4.999999866870943, 4.999999696844845, 3.8925884617955835, 4.999999846075734, 4.999999846079191, 3.6746483494128603, 2.5153283897091496, 4.346521214731182e-8, 4.346534512529267e-8, 4.3465342879947955e-8, 3.954083382104079e-8, 2.941606342979051e-8, 2.941602940532639e-8, 2.941602049179647e-8, 2.9415996475971363e-8, 2.941594466197134e-8, 2.9415844679636142e-8, 3.1201342908535565e-8, 3.12012000378666e-8, 3.9541480947856964e-8, 3.120129433516046e-8, 3.1201293194217276e-8, 3.120128891570256e-8, 4.666616176402652, 4.2965790094976954, 1.0288580673455925]), outer = 2)))	
	set_start_value.(b_plus, map(+, [0 for i in 1:length(T)+outage], 
		repeat(round.([0.39393066890710743, 5.990893352924856, 0.2075646242614788, 8.991289056801049e-7, 0.6201504854591727, 6.431780575016466, 8.990916082594305e-7, 8.991514072833661e-7, -6.3965450542933626e-9, 1.202109499351938e-6, 1.8292274104903435e-6, 2.6561597324247812e-6, 6.339605172051425, 9.045481795756041, 8.898024636588499, 7.935372820993773, 4.941131980929224, 0.5749666253542441, 8.99211013706141e-7, 8.991716948987293e-7, 1.8282767281342673e-7, 8.990899844163186e-7, 8.990897779606156e-7, 8.990905679292054e-7, 8.991017382754243e-7, 8.991060985822936e-7, 8.991163358188456e-7]), outer = 2)))	
	set_start_value.(b_minus, map(+, [0 for i in 1:length(T)+outage], repeat(round.([8.991829651168648e-7, 8.990437343271774e-7, 8.99094628390315e-7, 0.9607263363542307, 8.991747137671042e-7, 8.989750329840357e-7, 1.1833418871940131, 8.612234199778483, 8.502830005769852, 3.6264614166299318e-6, 1.7876113923673805e-6, 1.3694454383658838e-6, 8.991318835892509e-7, 8.991070221426838e-7, 8.991101774680498e-7, 8.991232399417761e-7, 8.991510286959842e-7, 8.992075728139885e-7, 0.2110259071714506, 4.83557189279799, 7.629696931387011, 6.847884500741222, 7.404882261596212, 5.190705198584314, 0.7692204722053042, 0.8943680694979523, 3.951231045319438]), outer = 2)))		
	set_start_value.(g_plus, map(+, [0 for i in 1:length(T)+outage], repeat(round.([0.38729850171589897, 5.908085022596629, 1.6524774435651045e-7, 8.152619482711186e-8, 0.25328955131771, 5.9247159468969945, 8.15254398225898e-8, 8.152667729623816e-8, 4.404924788957985, 8.489209514169394, 1.125312525608517, 6.650000850551692e-7, 6.649773711676688e-7, 6.649906727968432e-7, 6.64988972327515e-7, 6.649819638454448e-7, 6.649668035273976e-7, 6.649363873730011e-7, 2.214300287773459e-8, 2.2142959655748515e-8, 2.2142906964846808e-8, 2.214301516132642e-8, 2.3229968228832174e-8, 2.323001901552339e-8, 1.1757006791123573e-7, 9.146369931798415e-8, 7.423234367361756e-8]), outer = 2)))		
	set_start_value.(g_minus, map(+, [0 for i in 1:length(T)+outage], repeat(round.([1.4976942963122353e-8, 1.611214725097671e-8, 2.241559590697754e-8, 3.186030378027718e-8, 1.8723490087802263e-8, 1.8723289182528133e-8, 3.186046171422308e-8, 3.186020286242547e-8, 6.650212340168904e-7, 6.650000378372414e-7, 6.650003957256448e-7, 4.804144491931031, 4.0615886562564, 5.469226488530689, 7.905198209163032, 8.948508067394656, 10.670888875935118, 9.636747875130435, 5.54214611759292, 2.1677580901230695, 0.6766405508509762, 2.599888492323188e-7, 2.0179723873521621e-7, 2.017951758063469e-7, 2.4822911015552774e-8, 2.5161058185470468e-8, 2.5504917213530083e-8]), outer = 2)))	
	set_start_value.(r_Δ, map(+, [0 for i in 1:length(T)+outage], repeat(round.([1.2476226462335962e-7, 1.3395046558741293e-7, 3.096069868381261e-7, 1.1074115510221632, 1.5605666312287994e-7, 1.560522406585928e-7, 1.3253516634049065, 2.484671623108455, 5.000000003429293, 4.999999983590228, 4.999999964576983, 1.0418034683148722, 1.8285911781323915e-7, 1.8285721734624528e-7, 1.8285714277278304e-7, 1.8285713997651174e-7, 1.8285713721219306e-7, 1.8285722691803912e-7, 2.3473757729562566e-7, 2.3473504204850956e-7, 1.681894722765906, 2.347351513992217e-7, 1.859366421838519e-7, 1.8593477986377135e-7, 2.512451436708406e-7, 3.8190373068981424e-7, 7.735801223380891e-7]), outer = 2)))			
	set_start_value.(r_b_minus, map(+, [0 for i in 1:length(T)+outage], repeat(round.([1.4971710806072909, 7.360446212258328, 1.3170369204783252, 2.528740613223969e-6, 1.497169298678696, 7.360447786836681, 1.3170369759461422, 1.0582408189729051e-7, 1.497170090612351, 7.360449643643469, 1.3170370339620125, 1.4455210317281813e-7, 1.227023480750986e-5, 1.2262712573548463e-5, 1.2262434605322195e-5, 1.2262429118480203e-5, 1.2262386088277289e-5, 1.2262655572470277e-5, 1.7176578055847507, 4.114434767058682, 2.370302923876297, 1.9722605258651373, 1.6130359148677775, 1.4477464132079343, 1.0322141206823812, 1.1707738494132685, 1.3126957284486467]), outer = 2)))	
	set_start_value.(SOC, map(+, [0 for i in 1:length(T)+outage], repeat(round.([35.374233281237935, 41.065581112425, 41.26276665133351, 40.35007748596945, 40.939219592939686, 47.04941028517905, 45.92523634648176, 37.74361471088604, 29.665926199327966, 29.665923896193647, 29.665923935728863, 29.665925158107438, 35.688549217381, 44.28175606919757, 52.734878619801975, 60.273481945578986, 64.96755647326827, 65.51377391310761, 65.31330015554519, 60.719507711600215, 53.471295800468845, 46.96580637890017, 39.931169084519055, 35.0, 34.26924140555162, 33.419592593679354, 29.665923954786404]), outer = 2)))	

	### Define the objective function
	@objective(
		m, 
		Min, 
		(sum(inputs.c_Δ .* Δ
			+ inputs.c_b * (b_plus + b_minus)
			- inputs.p_d * inputs.d)
			+ (sum([sum([inputs.c_g_plus[t] * g_plus[t] - inputs.p_g_minus[t] * g_minus[t] + E[t] 
									   for t in 1:length(T)+inputs.κ if t ∉ τ:τ+inputs.κ])
									   for τ in 1:length(T)]) 
			    + (sum([sum([inputs.c_Δ .* r_Δ[t] + inputs.c_b .* r_b_minus[t] 
									   for t in τ:τ+inputs.κ]) 
									   for τ in eachindex(T)]))) * (inputs.ω  / (length(T)+inputs.κ)) 
			+ sum([inputs.c_g_plus[t] * g_plus[t] - inputs.p_g_minus[t] * g_minus[t] + E[t] 
									   for t in 1:length(T)+inputs.κ]) * (1-inputs.ω)
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
	for t in 1:length(T)+inputs.κ
		for τ in 1:max(t-inputs.κ, 1)
			@constraint(m, inputs.SOC_min <= SOC[t] - sum(inputs.η_minus * r_b_minus[s] * Δt / b_cap for s in τ:min(τ+inputs.κ, t))
						)
		end  
	end
	@constraint(m, [t = 1:length(T)+inputs.κ], 
				SOC[t] .- inputs.SOC_max .<= 0.00
				)
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

	x = r_b_minus + r_Δ + inputs.s_max + Δ + b_minus - b_plus - inputs.d

	for j in eachindex(T)
		register(m, Symbol("mvncdf_$j"), length(x), distribution(j, inputs)[1], distribution(j, inputs)[2])
		add_nonlinear_constraint(m, :($(Symbol("mvncdf_$j"))($(x...)) >= $(inputs.p)))
	end

	function define(t, inputs)
		c = inputs.c_g_plus[t] + c_γ_plus
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
	y = inputs.d - inputs.s_max - Δ - g_plus - b_minus + b_plus + g_minus
	for t in eachindex(y)
		register(m, Symbol("expected_$t"), length(y), define(t, inputs)[1], define(t, inputs)[2])
		add_nonlinear_constraint(m, :($(Symbol("expected_$t"))($(y...)) == $(E[t])))
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
	@variable(m, E[1:length(T)+inputs.κ] >= 0.)

	### add start values
	set_start_value.(Δ, map(+, [0 for i in 1:length(T)+outage], 
		repeat(round.([4.999999875405635, 4.999999866870943, 4.999999696844845, 3.8925884617955835, 4.999999846075734, 4.999999846079191, 3.6746483494128603, 2.5153283897091496, 4.346521214731182e-8, 4.346534512529267e-8, 4.3465342879947955e-8, 3.954083382104079e-8, 2.941606342979051e-8, 2.941602940532639e-8, 2.941602049179647e-8, 2.9415996475971363e-8, 2.941594466197134e-8, 2.9415844679636142e-8, 3.1201342908535565e-8, 3.12012000378666e-8, 3.9541480947856964e-8, 3.120129433516046e-8, 3.1201293194217276e-8, 3.120128891570256e-8, 4.666616176402652, 4.2965790094976954, 1.0288580673455925]), outer = 2)))	
	set_start_value.(b_plus, map(+, [0 for i in 1:length(T)+outage], 
		repeat(round.([0.39393066890710743, 5.990893352924856, 0.2075646242614788, 8.991289056801049e-7, 0.6201504854591727, 6.431780575016466, 8.990916082594305e-7, 8.991514072833661e-7, -6.3965450542933626e-9, 1.202109499351938e-6, 1.8292274104903435e-6, 2.6561597324247812e-6, 6.339605172051425, 9.045481795756041, 8.898024636588499, 7.935372820993773, 4.941131980929224, 0.5749666253542441, 8.99211013706141e-7, 8.991716948987293e-7, 1.8282767281342673e-7, 8.990899844163186e-7, 8.990897779606156e-7, 8.990905679292054e-7, 8.991017382754243e-7, 8.991060985822936e-7, 8.991163358188456e-7]), outer = 2)))	
	set_start_value.(b_minus, map(+, [0 for i in 1:length(T)+outage], repeat(round.([8.991829651168648e-7, 8.990437343271774e-7, 8.99094628390315e-7, 0.9607263363542307, 8.991747137671042e-7, 8.989750329840357e-7, 1.1833418871940131, 8.612234199778483, 8.502830005769852, 3.6264614166299318e-6, 1.7876113923673805e-6, 1.3694454383658838e-6, 8.991318835892509e-7, 8.991070221426838e-7, 8.991101774680498e-7, 8.991232399417761e-7, 8.991510286959842e-7, 8.992075728139885e-7, 0.2110259071714506, 4.83557189279799, 7.629696931387011, 6.847884500741222, 7.404882261596212, 5.190705198584314, 0.7692204722053042, 0.8943680694979523, 3.951231045319438]), outer = 2)))		
	set_start_value.(g_plus, map(+, [0 for i in 1:length(T)+outage], repeat(round.([0.38729850171589897, 5.908085022596629, 1.6524774435651045e-7, 8.152619482711186e-8, 0.25328955131771, 5.9247159468969945, 8.15254398225898e-8, 8.152667729623816e-8, 4.404924788957985, 8.489209514169394, 1.125312525608517, 6.650000850551692e-7, 6.649773711676688e-7, 6.649906727968432e-7, 6.64988972327515e-7, 6.649819638454448e-7, 6.649668035273976e-7, 6.649363873730011e-7, 2.214300287773459e-8, 2.2142959655748515e-8, 2.2142906964846808e-8, 2.214301516132642e-8, 2.3229968228832174e-8, 2.323001901552339e-8, 1.1757006791123573e-7, 9.146369931798415e-8, 7.423234367361756e-8]), outer = 2)))		
	set_start_value.(g_minus, map(+, [0 for i in 1:length(T)+outage], repeat(round.([1.4976942963122353e-8, 1.611214725097671e-8, 2.241559590697754e-8, 3.186030378027718e-8, 1.8723490087802263e-8, 1.8723289182528133e-8, 3.186046171422308e-8, 3.186020286242547e-8, 6.650212340168904e-7, 6.650000378372414e-7, 6.650003957256448e-7, 4.804144491931031, 4.0615886562564, 5.469226488530689, 7.905198209163032, 8.948508067394656, 10.670888875935118, 9.636747875130435, 5.54214611759292, 2.1677580901230695, 0.6766405508509762, 2.599888492323188e-7, 2.0179723873521621e-7, 2.017951758063469e-7, 2.4822911015552774e-8, 2.5161058185470468e-8, 2.5504917213530083e-8]), outer = 2)))	
	set_start_value.(r_Δ, map(+, [0 for i in 1:length(T)+outage], repeat(round.([1.2476226462335962e-7, 1.3395046558741293e-7, 3.096069868381261e-7, 1.1074115510221632, 1.5605666312287994e-7, 1.560522406585928e-7, 1.3253516634049065, 2.484671623108455, 5.000000003429293, 4.999999983590228, 4.999999964576983, 1.0418034683148722, 1.8285911781323915e-7, 1.8285721734624528e-7, 1.8285714277278304e-7, 1.8285713997651174e-7, 1.8285713721219306e-7, 1.8285722691803912e-7, 2.3473757729562566e-7, 2.3473504204850956e-7, 1.681894722765906, 2.347351513992217e-7, 1.859366421838519e-7, 1.8593477986377135e-7, 2.512451436708406e-7, 3.8190373068981424e-7, 7.735801223380891e-7]), outer = 2)))			
	set_start_value.(r_b_minus, map(+, [0 for i in 1:length(T)+outage], repeat(round.([1.4971710806072909, 7.360446212258328, 1.3170369204783252, 2.528740613223969e-6, 1.497169298678696, 7.360447786836681, 1.3170369759461422, 1.0582408189729051e-7, 1.497170090612351, 7.360449643643469, 1.3170370339620125, 1.4455210317281813e-7, 1.227023480750986e-5, 1.2262712573548463e-5, 1.2262434605322195e-5, 1.2262429118480203e-5, 1.2262386088277289e-5, 1.2262655572470277e-5, 1.7176578055847507, 4.114434767058682, 2.370302923876297, 1.9722605258651373, 1.6130359148677775, 1.4477464132079343, 1.0322141206823812, 1.1707738494132685, 1.3126957284486467]), outer = 2)))	
	set_start_value.(SOC, map(+, [0 for i in 1:length(T)+outage], repeat(round.([35.374233281237935, 41.065581112425, 41.26276665133351, 40.35007748596945, 40.939219592939686, 47.04941028517905, 45.92523634648176, 37.74361471088604, 29.665926199327966, 29.665923896193647, 29.665923935728863, 29.665925158107438, 35.688549217381, 44.28175606919757, 52.734878619801975, 60.273481945578986, 64.96755647326827, 65.51377391310761, 65.31330015554519, 60.719507711600215, 53.471295800468845, 46.96580637890017, 39.931169084519055, 35.0, 34.26924140555162, 33.419592593679354, 29.665923954786404]), outer = 2)))	

	### Define the objective function
	@objective(
		m, 
		Min, 
		(sum(inputs.c_Δ .* Δ
			+ inputs.c_b * (b_plus + b_minus)
			- inputs.p_d * inputs.d)
			+ (sum([sum([inputs.c_g_plus[t] * g_plus[t] - inputs.p_g_minus[t] * g_minus[t] + E[t] 
									   for t in 1:length(T)+inputs.κ if t ∉ τ:τ+inputs.κ])
									   for τ in 1:length(T)]) 
			    + (sum([sum([inputs.c_Δ .* r_Δ[t] + inputs.c_b .* r_b_minus[t] 
									   for t in τ:τ+inputs.κ]) 
									   for τ in eachindex(T)]))) * (inputs.ω  / (length(T)+inputs.κ)) 
			+ sum([inputs.c_g_plus[t] * g_plus[t] - inputs.p_g_minus[t] * g_minus[t] + E[t] 
									   for t in 1:length(T)+inputs.κ]) * (1-inputs.ω)
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
	for t in 1:length(T)+inputs.κ
		for τ in 1:max(t-inputs.κ, 1)
			@constraint(m, inputs.SOC_min <= SOC[t] - sum(inputs.η_minus * r_b_minus[s] * Δt / b_cap for s in τ:min(τ+inputs.κ, t))
						)
		end  
	end
	@constraint(m, [t = 1:length(T)+inputs.κ], 
				SOC[t] .- inputs.SOC_max .<= 0.00
				)
	@constraint(m, SOC[length(T)] == inputs.SOC_0)

	function create_cons(w, inputs)
		k = (1+inputs.κ) * (w-1) + 1
		f(x...) = ProbFunction(vec([x[r] for r in w:w+inputs.κ]), SampleOnSphere[:, k:k+inputs.κ], inputs.μ[w:w+inputs.κ])[1]
		function ∇f(g::AbstractVector{T}, x::T...) where {T}
			for i in eachindex(x)
				if i in w:w+inputs.κ
					g[i] = ProbFunction(vec([x[r] for r in w:w+inputs.κ]), SampleOnSphere[:, k:k+inputs.κ], inputs.μ[w:w+inputs.κ])[2][i-w+1]
				else 
					g[i] = 0.00000
				end 
			end
			return
		end
		return f, ∇f
	end 
	
	x = r_b_minus + r_Δ + inputs.s_max + Δ + b_minus - b_plus - inputs.d

	for j in eachindex(T)
		register(m, Symbol("srd_prob_$j"), length(x), create_cons(j, inputs)[1], create_cons(j, inputs)[2])
		add_nonlinear_constraint(m, :($(Symbol("srd_prob_$j"))($(x...)) >= $(inputs.p)))
	end

	function define(t, inputs)
		# norm_pdf(k) = exp(-(k^2) / 2) / (2 * pi)^0.5
		c = inputs.c_g_plus[t] + c_γ_plus
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
	y = inputs.d - inputs.s_max - Δ - g_plus - b_minus + b_plus + g_minus
	for t in eachindex(y)
		register(m, Symbol("expected_$t"), length(y), define(t, inputs)[1], define(t, inputs)[2])
		add_nonlinear_constraint(m, :($(Symbol("expected_$t"))($(y...)) == $(E[t])))
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
	@variable(m, 0. <= p <= 1., start = 0.9)
	@variable(m, E[1:length(T)+inputs.κ] >= 0., start = 0.00)

	### add start values
	set_start_value.(Δ, map(+, [0 for i in 1:length(T)+outage], 
		repeat(round.([4.999999875405635, 4.999999866870943, 4.999999696844845, 3.8925884617955835, 4.999999846075734, 4.999999846079191, 3.6746483494128603, 2.5153283897091496, 4.346521214731182e-8, 4.346534512529267e-8, 4.3465342879947955e-8, 3.954083382104079e-8, 2.941606342979051e-8, 2.941602940532639e-8, 2.941602049179647e-8, 2.9415996475971363e-8, 2.941594466197134e-8, 2.9415844679636142e-8, 3.1201342908535565e-8, 3.12012000378666e-8, 3.9541480947856964e-8, 3.120129433516046e-8, 3.1201293194217276e-8, 3.120128891570256e-8, 4.666616176402652, 4.2965790094976954, 1.0288580673455925]), outer = 2)))	
	set_start_value.(b_plus, map(+, [0 for i in 1:length(T)+outage], 
		repeat(round.([0.39393066890710743, 5.990893352924856, 0.2075646242614788, 8.991289056801049e-7, 0.6201504854591727, 6.431780575016466, 8.990916082594305e-7, 8.991514072833661e-7, -6.3965450542933626e-9, 1.202109499351938e-6, 1.8292274104903435e-6, 2.6561597324247812e-6, 6.339605172051425, 9.045481795756041, 8.898024636588499, 7.935372820993773, 4.941131980929224, 0.5749666253542441, 8.99211013706141e-7, 8.991716948987293e-7, 1.8282767281342673e-7, 8.990899844163186e-7, 8.990897779606156e-7, 8.990905679292054e-7, 8.991017382754243e-7, 8.991060985822936e-7, 8.991163358188456e-7]), outer = 2)))	
	set_start_value.(b_minus, map(+, [0 for i in 1:length(T)+outage], repeat(round.([8.991829651168648e-7, 8.990437343271774e-7, 8.99094628390315e-7, 0.9607263363542307, 8.991747137671042e-7, 8.989750329840357e-7, 1.1833418871940131, 8.612234199778483, 8.502830005769852, 3.6264614166299318e-6, 1.7876113923673805e-6, 1.3694454383658838e-6, 8.991318835892509e-7, 8.991070221426838e-7, 8.991101774680498e-7, 8.991232399417761e-7, 8.991510286959842e-7, 8.992075728139885e-7, 0.2110259071714506, 4.83557189279799, 7.629696931387011, 6.847884500741222, 7.404882261596212, 5.190705198584314, 0.7692204722053042, 0.8943680694979523, 3.951231045319438]), outer = 2)))		
	set_start_value.(g_plus, map(+, [0 for i in 1:length(T)+outage], repeat(round.([0.38729850171589897, 5.908085022596629, 1.6524774435651045e-7, 8.152619482711186e-8, 0.25328955131771, 5.9247159468969945, 8.15254398225898e-8, 8.152667729623816e-8, 4.404924788957985, 8.489209514169394, 1.125312525608517, 6.650000850551692e-7, 6.649773711676688e-7, 6.649906727968432e-7, 6.64988972327515e-7, 6.649819638454448e-7, 6.649668035273976e-7, 6.649363873730011e-7, 2.214300287773459e-8, 2.2142959655748515e-8, 2.2142906964846808e-8, 2.214301516132642e-8, 2.3229968228832174e-8, 2.323001901552339e-8, 1.1757006791123573e-7, 9.146369931798415e-8, 7.423234367361756e-8]), outer = 2)))		
	set_start_value.(g_minus, map(+, [0 for i in 1:length(T)+outage], repeat(round.([1.4976942963122353e-8, 1.611214725097671e-8, 2.241559590697754e-8, 3.186030378027718e-8, 1.8723490087802263e-8, 1.8723289182528133e-8, 3.186046171422308e-8, 3.186020286242547e-8, 6.650212340168904e-7, 6.650000378372414e-7, 6.650003957256448e-7, 4.804144491931031, 4.0615886562564, 5.469226488530689, 7.905198209163032, 8.948508067394656, 10.670888875935118, 9.636747875130435, 5.54214611759292, 2.1677580901230695, 0.6766405508509762, 2.599888492323188e-7, 2.0179723873521621e-7, 2.017951758063469e-7, 2.4822911015552774e-8, 2.5161058185470468e-8, 2.5504917213530083e-8]), outer = 2)))	
	set_start_value.(r_Δ, map(+, [0 for i in 1:length(T)+outage], repeat(round.([1.2476226462335962e-7, 1.3395046558741293e-7, 3.096069868381261e-7, 1.1074115510221632, 1.5605666312287994e-7, 1.560522406585928e-7, 1.3253516634049065, 2.484671623108455, 5.000000003429293, 4.999999983590228, 4.999999964576983, 1.0418034683148722, 1.8285911781323915e-7, 1.8285721734624528e-7, 1.8285714277278304e-7, 1.8285713997651174e-7, 1.8285713721219306e-7, 1.8285722691803912e-7, 2.3473757729562566e-7, 2.3473504204850956e-7, 1.681894722765906, 2.347351513992217e-7, 1.859366421838519e-7, 1.8593477986377135e-7, 2.512451436708406e-7, 3.8190373068981424e-7, 7.735801223380891e-7]), outer = 2)))			
	set_start_value.(r_b_minus, map(+, [0 for i in 1:length(T)+outage], repeat(round.([1.4971710806072909, 7.360446212258328, 1.3170369204783252, 2.528740613223969e-6, 1.497169298678696, 7.360447786836681, 1.3170369759461422, 1.0582408189729051e-7, 1.497170090612351, 7.360449643643469, 1.3170370339620125, 1.4455210317281813e-7, 1.227023480750986e-5, 1.2262712573548463e-5, 1.2262434605322195e-5, 1.2262429118480203e-5, 1.2262386088277289e-5, 1.2262655572470277e-5, 1.7176578055847507, 4.114434767058682, 2.370302923876297, 1.9722605258651373, 1.6130359148677775, 1.4477464132079343, 1.0322141206823812, 1.1707738494132685, 1.3126957284486467]), outer = 2)))	
	set_start_value.(SOC, map(+, [0 for i in 1:length(T)+outage], repeat(round.([35.374233281237935, 41.065581112425, 41.26276665133351, 40.35007748596945, 40.939219592939686, 47.04941028517905, 45.92523634648176, 37.74361471088604, 29.665926199327966, 29.665923896193647, 29.665923935728863, 29.665925158107438, 35.688549217381, 44.28175606919757, 52.734878619801975, 60.273481945578986, 64.96755647326827, 65.51377391310761, 65.31330015554519, 60.719507711600215, 53.471295800468845, 46.96580637890017, 39.931169084519055, 35.0, 34.26924140555162, 33.419592593679354, 29.665923954786404]), outer = 2)))	

	### Define the objective function
	@objective(
		m, 
		Max,
		p)

	### Define the constraints
	@constraint(m, Δ + r_Δ .- inputs.Δ_max .<= 0.00)
	@constraint(m, b_minus + r_b_minus .- inputs.b_max .<= 0.00)
	@constraint(m, [t = 2:length(T)+inputs.κ], 
					SOC[t] == SOC[t-1] + (inputs.η_plus * b_plus[t] - inputs.η_minus * b_minus[t]) * Δt / b_cap
				)
	@constraint(m, 
				SOC[1] == inputs.SOC_0 + (inputs.η_plus * b_plus[1] - inputs.η_minus * b_minus[1]) * Δt / b_cap
				)
	for t in 1:length(T)+inputs.κ
		for τ in 1:max(t-inputs.κ, 1)
			@constraint(m, inputs.SOC_min <= SOC[t] - sum(inputs.η_minus * r_b_minus[s] * Δt / b_cap for s in τ:min(τ+inputs.κ, t))
						)
		end  
	end
	@constraint(m, [t = 1:length(T)+inputs.κ], 
				SOC[t] .- inputs.SOC_max .<= 0.00
				)
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

	x = r_b_minus + r_Δ + inputs.s_max + Δ + b_minus - b_plus - inputs.d

	for j in eachindex(T)
		register(m, Symbol("mvncdf_$j"), length(x), distribution(j, inputs)[1], distribution(j, inputs)[2])
		add_nonlinear_constraint(m, :($(Symbol("mvncdf_$j"))($(x...)) >= $(p)))
	end

	function define(t, inputs)
		# norm_pdf(k) = exp(-(k^2) / 2) / (2 * pi)^0.5
		c = inputs.c_g_plus[t] + c_γ_plus
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
	y = inputs.d - inputs.s_max - Δ - g_plus - b_minus + b_plus + g_minus
	for t in eachindex(y)
		register(m, Symbol("expected_$t"), length(y), define(t, inputs)[1], define(t, inputs)[2])
		add_nonlinear_constraint(m, :($(Symbol("expected_$t"))($(y...)) == $(E[t])))
	end

	return m
end