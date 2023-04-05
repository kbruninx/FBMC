using JuMP, BilevelJuMP, Gurobi, Dualization
using DataFrames, XLSX
using LinearAlgebra
using Alpine
using Ipopt
using Statistics
using QuadraticToBinary
using Plots
using SparseArrays
using Formatting

function sum_z_np(np, num_t)
    result = []
    for t in 1:num_t
        sum = 0
        for z in 1:num_z
            sum += np[num_t*(z-1)+t]
        end
        push!(result, sum)
    end
    return result
end

fbmc_zones =  ["ALBE", "ALDE", "AT", "BE", "CZ", "DE_LU", "FR", "HR", "HU", "NL", "PL", "RO", "SI", "SK"]
non_fbmc_zones =  ["CH", "GB", "ES", "IT_NORD", "DK_1", "NO_2"]
all_zones = vcat(fbmc_zones, non_fbmc_zones)

function has_interconnector(z1, z2)
    if sort([all_zones[z1], all_zones[z2]]) == ["FR", "GB"]
        return true
    elseif sort([all_zones[z1], all_zones[z2]]) == ["ES", "FR"]
        return true
    elseif sort([all_zones[z1], all_zones[z2]]) == ["CH", "FR"]
        return true
    elseif sort([all_zones[z1], all_zones[z2]]) == ["FR", "IT_NORD"]
        return true
    elseif sort([all_zones[z1], all_zones[z2]]) == ["AT", "IT_NORD"]
        return true
    elseif sort([all_zones[z1], all_zones[z2]]) == ["AT", "CH"]
        return true
    elseif sort([all_zones[z1], all_zones[z2]]) == ["BE", "GB"]
        return true
    elseif sort([all_zones[z1], all_zones[z2]]) == ["DK_1", "NL"]
        return true
    elseif sort([all_zones[z1], all_zones[z2]]) == ["GB", "NL"]
        return true
    elseif sort([all_zones[z1], all_zones[z2]]) == ["NL", "NO_2"]
        return true
    elseif sort([all_zones[z1], all_zones[z2]]) == ["IT_NORD", "SI"]
        return true
end

function atc(z1, z2, t)

end

num_t_passed = 0
experiments = [120]
#experiments = [240, 240, 240, 240, 240, 240, 240, 240, 240, 240, 253, 253]

experiment_results_alpha = []
experiment_results_beta = []
experiment_results_gamma = []

# iterations out: 4, 6, 9, 10, 11
# infeasible: 6
"""
num_t_passed = 4*240
experiments = [240]
"""

alpha_global = zeros((num_z+num_z_non_fbmc)*num_tech)
beta_global = zeros((num_z+num_z_non_fbmc)*num_tech)
gamma_global = zeros((num_z+num_z_non_fbmc)*num_tech)

iteration_count = 1
for exp_len in experiments
    println("ITERATION ", iteration_count)
    
    num_t = exp_len

    coal_prices = coal_prices_g[(num_t_passed+1):(num_t_passed)+num_t]
    oil_prices = oil_prices_g[(num_t_passed+1):(num_t_passed)+num_t]
    gas_prices = gas_prices_g[(num_t_passed+1):(num_t_passed)+num_t]
    eua_prices = eua_prices_g[(num_t_passed+1):(num_t_passed)+num_t]

    at_obs = vec(at_obs_g[(num_t_passed+1):(num_t_passed)+num_t, :])
    be_obs = vec(be_obs_g[(num_t_passed+1):(num_t_passed)+num_t, :])
    cz_obs = vec(cz_obs_g[(num_t_passed+1):(num_t_passed)+num_t, :])
    de_obs = vec(de_obs_g[(num_t_passed+1):(num_t_passed)+num_t, :])
    fr_obs = vec(fr_obs_g[(num_t_passed+1):(num_t_passed)+num_t, :])
    hr_obs = vec(hr_obs_g[(num_t_passed+1):(num_t_passed)+num_t, :])
    hu_obs = vec(hu_obs_g[(num_t_passed+1):(num_t_passed)+num_t, :])
    nl_obs = vec(nl_obs_g[(num_t_passed+1):(num_t_passed)+num_t, :])
    pl_obs = vec(pl_obs_g[(num_t_passed+1):(num_t_passed)+num_t, :])
    ro_obs = vec(ro_obs_g[(num_t_passed+1):(num_t_passed)+num_t, :])
    si_obs = vec(si_obs_g[(num_t_passed+1):(num_t_passed)+num_t, :])
    sk_obs = vec(sk_obs_g[(num_t_passed+1):(num_t_passed)+num_t, :])
    albe_obs = zeros(num_t*num_tech)
    alde_obs = zeros(num_t*num_tech)

    ch_obs = vec(ch_obs_g[(num_t_passed+1):(num_t_passed)+num_t, :])
    gb_obs = vec(gb_obs_g[(num_t_passed+1):(num_t_passed)+num_t, :])
    es_obs = vec(es_obs_g[(num_t_passed+1):(num_t_passed)+num_t, :])
    it_obs = vec(it_obs_g[(num_t_passed+1):(num_t_passed)+num_t, :])
    dk_obs = vec(dk_obs_g[(num_t_passed+1):(num_t_passed)+num_t, :])
    no_obs = vec(no_obs_g[(num_t_passed+1):(num_t_passed)+num_t, :])

    g_obs = vcat(albe_obs, alde_obs, at_obs, be_obs, cz_obs, de_obs, fr_obs, hr_obs, hu_obs, nl_obs, pl_obs, ro_obs, si_obs, sk_obs, ch_obs, gb_obs, es_obs, it_obs, dk_obs, no_obs)
    g_obs = convert(Vector{Float64}, g_obs)

    demand_fbmc = vec(demand_g[(num_t_passed+1):(num_t_passed)+num_t, :])
    demand_non_fbmc = vec(demand_g_non_fbmc[(num_t_passed+1):(num_t_passed)+num_t, :])
    demand = vcat(demand_fbmc, demand_non_fbmc)

    lambda_obs_fbmc = vec(lambda_obs_g[(num_t_passed+1):(num_t_passed)+num_t, :])
    lambda_obs_non_fbmc = vec(lambda_obs_g_non_fbmc[(num_t_passed+1):(num_t_passed)+num_t, :])
    lambda_obs = vcat(lambda_obs_fbmc, lambda_obs_non_fbmc)

    np_obs = vec(np_obs_g[(num_t_passed+1):(num_t_passed)+num_t, :])

    demand = convert(Vector{Float64}, demand)
    lambda_obs = convert(Vector{Float64}, lambda_obs)
    np_obs = convert(Vector{Float64}, np_obs)

    ptdf_z = vec(ptdf_z_g[(num_t_passed*num_j+1):(num_t_passed*num_j)+num_t*num_j, :])
    ram = vec(ram_g[(num_t_passed*num_j+1):(num_t_passed*num_j)+num_t*num_j, :])
    ptdf_z = convert(Vector{Float64}, ptdf_z)
    ram = convert(Vector{Float64}, ram)

    ren_gen_fbmc = vec(ren_gen_g[(num_t_passed+1):(num_t_passed)+num_t, :])
    ren_gen_non_fbmc = vec(ren_gen_g_non_fbmc[(num_t_passed+1):(num_t_passed)+num_t, :])
    ren_gen = vcat(ren_gen_fbmc, ren_gen_non_fbmc)

    ren_gen = convert(Vector{Float64}, ren_gen)

    g_max_t_fbmc = zeros(num_z*num_tech*num_t)
    for z in 1:num_z
        for t in 1:num_t
            for tech in 1:num_tech
                g_max_t_fbmc[num_t*num_tech*(z-1)+num_t*(tech-1)+t] = g_max_g[num_tech*(z-1)+tech] 
            end
        end
    end

    g_max_t_non_fbmc = zeros(num_z_non_fbmc*num_tech*num_t)
    for z in 1:num_z_non_fbmc
        for t in 1:num_t
            for tech in 1:num_tech
                g_max_t_non_fbmc[num_t*num_tech*(z-1)+num_t*(tech-1)+t] = g_max_g_non_fbmc[num_tech*(z-1)+tech] 
            end
        end
    end

    g_max_t = vcat(g_max_t_fbmc, g_max_t_non_fbmc)

    model = direct_model(Gurobi.Optimizer())
    set_time_limit_sec(model, 120.0)
    set_optimizer_attribute(model, "NonConvex", 2)

    @variable(model, c[1:(num_z+num_z_non_fbmc)*num_tech*num_t] >= 0)

    @variable(model, alpha[1:(num_z+num_z_non_fbmc)*num_tech])
    @variable(model, beta[1:(num_z+num_z_non_fbmc)*num_tech])
    @variable(model, gamma[1:(num_z+num_z_non_fbmc)*num_tech])

    @constraint(model, alpha .+ alpha_global .>= 0)
    @constraint(model, beta .+ beta_global .>= 0)
    @constraint(model, gamma .+ gamma_global .>= 0)

    for z in 1:(num_z+num_z_non_fbmc)
        for tech in 3:num_tech
            if tech == 7
                @constraint(model, gamma_global[num_tech*(z-1)+tech] + gamma[num_tech*(z-1)+tech] == 0)
            end
        end
        for t in 1:num_t
            for tech in 1:num_tech

                alpha_coeff = alpha_global .+ alpha
                beta_coeff = beta_global .+ beta
                gamma_coeff = gamma_global .+ gamma

                if tech == 1 || tech == 8 || tech == 9 || tech == 10 # without fuel-responsiveness
                    @constraint(model, c[num_t*num_tech*(z-1)+num_t*(tech-1)+t] == alpha_coeff[num_tech*(z-1)+tech] + (beta_coeff[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t] + gamma_coeff[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t]^2))
                elseif tech == 2 || tech == 3 || tech == 5 # coal
                    @constraint(model, c[num_t*num_tech*(z-1)+num_t*(tech-1)+t] == coal_prices[t]/10 * (alpha_coeff[num_tech*(z-1)+tech] + beta_coeff[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t] + gamma_coeff[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t]^2))
                elseif tech == 4 # gas
                    @constraint(model, c[num_t*num_tech*(z-1)+num_t*(tech-1)+t] == gas_prices[t]/10 * (alpha_coeff[num_tech*(z-1)+tech] + beta_coeff[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t] + gamma_coeff[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t]^2))
                elseif tech == 6 # oil
                    @constraint(model, c[num_t*num_tech*(z-1)+num_t*(tech-1)+t] == oil_prices[t]/10 * (alpha_coeff[num_tech*(z-1)+tech] + beta_coeff[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t] + gamma_coeff[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t]^2))
                elseif tech == 7 # hydro
                    @constraint(model, c[num_t*num_tech*(z-1)+num_t*(tech-1)+t] == alpha_coeff[num_tech*(z-1)+tech] + beta_coeff[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t])
                end
            end
        end
    end

    @variable(model, eps1[1:(num_z+num_z_non_fbmc)*num_t] >= 0)
    @variable(model, eps2[1:(num_z+num_z_non_fbmc)*num_t] >= 0)
    #@variable(model, eps3[1:num_z*num_t] >= 0)
    #@variable(model, eps4[1:num_z*num_t] >= 0)

    @variable(model, g[1:(num_z+num_z_non_fbmc)*num_tech*num_t] >= 0)
    @variable(model, np[1:num_z*num_t])

    # dual variables
    @variable(model, lambda[1:(num_z+num_z_non_fbmc)*num_t])
    @variable(model, lambda_exchange[1:num_t])

    @variable(model, mu_gen[1:(num_z+num_z_non_fbmc)*num_tech*num_t] <= 0)
    @variable(model, mu_exchange[1:(num_j*num_t+2*(num_z+num_z_non_fbmc)*num_t)] <= 0)

    @constraint(model, g .== g_obs) # observe generation

    A_balance = spzeros((num_z+num_z_non_fbmc)*num_t, (num_z+num_z_non_fbmc)*num_tech*num_t+num_z*num_t+2*(num_z+num_z_non_fbmc)*num_t) # contains g and np and atc (e/i)
    prev_pos = (num_z+num_z_non_fbmc)*num_tech*num_t
    for z in 1:(num_z+num_z_non_fbmc)
        for t in 1:num_t
            for tech in 1:num_tech
                A_balance[num_t*(z-1)+t, num_t*num_tech*(z-1)+num_t*(tech-1)+t] = 1
            end
        end
    end

    for z in 1:num_z
        for t in 1:num_t
            A_balance[num_t*(z-1)+t, prev_pos+num_t*(z-1)+t] = -1 # np
        end
    end

    prev_pos += num_z*num_t

    for z1 in 1:(num_z+num_z_non_fbmc)
        for z2 in 1:(num_z+num_z_non_fbmc)
            for t in 1:num_t
                A_balance[num_t*(z1-1)+t, prev_pos+num_t*(z2-1)+t] = -1*has_interconnector(z1, z2) # exports
            end
        end
    end

    prev_pos += (num_z+num_z_non_fbmc)*num_t

    for z1 in 1:(num_z+num_z_non_fbmc)
        for z2 in 1:(num_z+num_z_non_fbmc)
            for t in 1:num_t
                A_balance[num_t*(z1-1)+t, prev_pos+num_t*(z2-1)+t] = has_interconnector(z1, z2) # imports
            end
        end
    end

    b1_balance = demand - ren_gen

    B_gen = sparse(cat(Matrix(I, (num_z+num_z_non_fbmc)*num_tech*num_t, (num_z+num_z_non_fbmc)*num_tech*num_t), spzeros((num_z+num_z_non_fbmc)*num_tech*num_t, num_z*num_t); dims=(2)))
    b2_gen = g_max_t

    A_exchange = spzeros(num_t, num_z*num_tech*num_t+num_z*num_t+2*(num_z+num_z_non_fbmc)*num_t)

    prev_pos = (num_z+num_z_non_fbmc)*num_tech*num_t
    for t in 1:num_t
        for z in 1:num_z
            A_exchange[t, prev_pos+num_t*(z-1)+t] = 1
        end
    end

    b1_exchange = spzeros(num_t)

    B_exchange_temp = spzeros(num_j*num_t, num_z*num_t)
    b2_exchange = ram
    for t in 1:num_t
        for j in 1:num_j
            for z in 1:num_z
                B_exchange_temp[num_t*(j-1)+t, num_t*(z-1)+t] = ptdf_z[num_j*num_t*(z-1) + num_j*(t-1) + j]
            end
        end
    end
    B_exchange_ptdf = sparse(cat(spzeros(num_j*num_t, (num_z+num_z_non_fbmc)*num_tech*num_t), B_exchange_temp; dims=(2)))
    
    B_exchange_atc = spzeros((num_z+num_z_non_fbmc)*num_t, (num_z+num_z_non_fbmc)*num_tech*num_t + num_z*num_t + 2*(num_z+num_z_non_fbmc)*num_t)
    prev_pos = (num_z+num_z_non_fbmc)*num_tech*num_t + num_z*num_t
    
    for z1 in 1:(num_z+num_z_non_fbmc)
        for z2 in 1:(num_z+num_z_non_fbmc)
            for t in 1:num_t
                B_exchange_atc[num_t*(z1-1)+t, prev_pos+num_t*(z2-1)+t] = atc(z1, z2, t) # exports
            end
        end
    end

    prev_pos += (num_z+num_z_non_fbmc)*num_t

    for z1 in 1:(num_z+num_z_non_fbmc)
        for z2 in 1:(num_z+num_z_non_fbmc)
            for t in 1:num_t
                B_exchange_atc[num_t*(z1-1)+t, prev_pos+num_t*(z2-1)+t] = atc(z2, z1, t) # imports
            end
        end
    end

    B_exchange = sparse(cat(B_exchange_ptdf, B_exchange_atc; dims=(1)))

    # primal constraints
    #@constraint(model, balance, A_balance * vcat(g, np) .== b1_balance) # demand equality
    @constraint(model, B_exchange * vcat(g, np) .<= b2_exchange) # net position, rams and ptdfs
    @constraint(model, sum_z_np(np, num_t) .== 0)

    # dual constraints
    @constraint(model, cat(A_balance, A_exchange; dims=(1))' * vcat(lambda, lambda_exchange) .+ cat(B_gen, B_exchange; dims=(1))' * vcat(mu_gen, mu_exchange) .== vcat(c, spzeros(num_z*num_t)))

    # strong duality gap theorem
    @constraint(model, sum(c) .- b1_balance' * lambda .- b1_exchange' * lambda_exchange .- b2_gen' * mu_gen .- b2_exchange' * mu_exchange == 0)

    @constraint(model, lambda .- lambda_obs .== eps1 .- eps2)
    #@constraint(model, np .- np_obs .== eps3 .- eps4)

    u = ones(num_z*num_t)
    @objective(model, Min, eps1' * u + eps2' * u)
    #@objective(model, Min, (eps1' * u + eps2' * u)/(num_z*num_t*maximum(lambda_obs)) + (eps3' * u + eps4' * u)/(num_z*num_t*maximum(np_obs)))

    optimize!(model)

    # iteratively adjust upon the previous values
    global num_t_passed += exp_len 
    #global alpha_global += JuMP.value.(alpha)
    #global beta_global += JuMP.value.(beta)
    #global gamma_global += JuMP.value.(gamma)

    if termination_status(model) == OPTIMAL
        println("ITERATION ", iteration_count, " FOUND OPTIMAL SOLUTION")
        push!(experiment_results_alpha, JuMP.value.(alpha))
        push!(experiment_results_beta, JuMP.value.(beta))
        push!(experiment_results_gamma, JuMP.value.(gamma))
    end
   
    global iteration_count += 1  
end

# SAVE COEFFICIENTS

df_coeffs = []
for z in 3:num_z
    alpha_coeffs = zeros(num_tech)
    beta_coeffs = zeros(num_tech)
    gamma_coeffs = zeros(num_tech)

    for tech in 1:num_tech
        alpha_exp = []
        beta_exp = []
        gamma_exp = []

        for e in 1:size(experiment_results_alpha)[1]
            push!(alpha_exp, experiment_results_alpha[e][num_tech*(z-1)+tech])
            push!(beta_exp, experiment_results_beta[e][num_tech*(z-1)+tech])
            push!(gamma_exp, experiment_results_gamma[e][num_tech*(z-1)+tech])
        end

        alpha_exp = convert(Vector{Float64}, alpha_exp)
        beta_exp = convert(Vector{Float64}, beta_exp)
        gamma_exp = convert(Vector{Float64}, gamma_exp)

        alpha_mean = mean(filter(!iszero, alpha_exp))
        beta_mean = mean(filter(!iszero, beta_exp))
        gamma_mean = mean(filter(!iszero, gamma_exp))

        if !isnan(alpha_mean)
            alpha_coeffs[tech] = alpha_mean
        end
        if !isnan(beta_mean)
            beta_coeffs[tech] = beta_mean
        end
        if !isnan(gamma_mean)
            gamma_coeffs[tech] = gamma_mean
        end
    end
    
    push!(df_coeffs, DataFrames.DataFrame(alpha=alpha_coeffs, beta=beta_coeffs, gamma=gamma_coeffs))
end

XLSX.writetable("cost_coefficients_alpha_fuel.xlsx",
    "AT" => df_coeffs[1],
    "BE" => df_coeffs[2],
    "CZ" => df_coeffs[3],
    "DE_LU" => df_coeffs[4],
    "FR" => df_coeffs[5],
    "HR" => df_coeffs[6],
    "HU" => df_coeffs[7],
    "NL" => df_coeffs[8],
    "PL" => df_coeffs[9],
    "RO" => df_coeffs[10],
    "SI" => df_coeffs[11],
    "SK" => df_coeffs[12],
    overwrite=true
)

"""
df_zones = []
for z in 3:num_z
    date_time = []
    plant_type = []
    marginal_cost = []
    production_level = []

    i = 1
    for t in 1:num_t
        for tech in 1:num_tech
            if JuMP.value.(c)[num_t*num_tech*(z-1)+num_t*(tech-1)+t] > 0 && JuMP.value.(g)[num_t*num_tech*(z-1)+num_t*(tech-1)+t] > 0
                push!(date_time, t)
                push!(plant_type, plant_types[tech])
                push!(marginal_cost, JuMP.value.(c)[num_t*num_tech*(z-1)+num_t*(tech-1)+t])
                push!(production_level, JuMP.value.(g)[num_t*num_tech*(z-1)+num_t*(tech-1)+t])
                i += 1
            end
        end
    end

    push!(df_zones, DataFrames.DataFrame(date_time=date_time, plant_type=plant_type, marginal_cost=marginal_cost, production_level=production_level))
end

XLSX.writetable("marginal_prices_output.xlsx",
    "AT" => df_zones[1],
    "BE" => df_zones[2],
    "CZ" => df_zones[3],
    "DE_LU" => df_zones[4],
    "FR" => df_zones[5],
    "HR" => df_zones[6],
    "HU" => df_zones[7],
    "NL" => df_zones[8],
    "PL" => df_zones[9],
    "RO" => df_zones[10],
    "SI" => df_zones[11],
    "SK" => df_zones[12],
    overwrite=true
)
"""