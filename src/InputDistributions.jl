## File for reading the input distributions generated by the ARIMA model

if sims == 10
    d = XLSX.readdata("$dir/CovMat/Data_10Simulations/D_forecast10sims.xlsx", "D_forecast10sims_two", "A2:M$(lpad(length(T)+outage+1, length(digits(Int(length(T)+outage+1)))))")
    s_dry = XLSX.readdata("$dir/CovMat/Data_10Simulations/S_dry_forecast10.xlsx", "S_dry_forecast10", "A2:M$(lpad(length(T)+outage+1, length(digits(Int(length(T)+outage+1)))))")

    d_error = transpose(XLSX.readdata("$dir/CovMat/Data_10Simulations/D_delta10sims.xlsx", "D_delta10sims_two", "A2:X11")) # TODO update this
    s_dry_error = XLSX.readdata("$dir/CovMat/Data_10Simulations/S_dry_error_xi10.xlsx", "S_dry_error_xi10", "A2:J$(lpad(length(T)+outage+1, length(digits(Int(length(T)+outage+1)))))")
    
else
    ## 100 simulations
    d = XLSX.readdata("$dir/CovMat/Data_100Simulations/D_forecast100sims_two.xlsx", "D_forecast100sims_two", "A2:CY$(lpad(length(T)+outage+1, length(digits(Int(length(T)+outage+1)))))")
    s_dry = XLSX.readdata("$dir/CovMat/Data_100Simulations/S_dry_forecast100.xlsx", "S_dry_forecast100", "A2:CY$(lpad(length(T)+outage+1, length(digits(Int(length(T)+outage+1)))))")

    d_error = transpose(XLSX.readdata("$dir/CovMat/Data_100Simulations/D_delta100sims_two.xlsx", "D_delta100sims_two", "A2:AY101")) # TODO generalize AA instead of selection in line 37 for d_error
    s_dry_error = XLSX.readdata("$dir/CovMat/Data_100Simulations/S_dry_error_xi100.xlsx", "S_dry_error_xi100", "A2:CV$(lpad(length(T)+outage+1, length(digits(Int(length(T)+outage+1)))))")
end

load_mean = [mean(d_error[:,i]) for i in 1:length(T)+outage]
drySun_mean = [mean(s_dry_error[:,i]) for i in 1:length(T)+outage]

# errors are identically and independently distributed:
CovMatrix = cov(transpose(s_dry_error)) + cov(transpose(d_error[1:length(T)+outage,:]))
Mean =  load_mean + drySun_mean

μ = vec([0 for i in 1:length(T)+outage])
Σ = CovMatrix
σ = diag(Σ, 0).^0.5

# Error function for the individual case (defined with the variance)
ξ = []
for i in 1:length(T)+outage
	append!(ξ, [Normal(μ[i], σ[i])])
end

a = vec([-Inf for i in 1:length(T)+outage])

