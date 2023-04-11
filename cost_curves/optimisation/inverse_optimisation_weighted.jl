using JuMP, HiGHS
#using Gurobi
using DataFrames, XLSX
using LinearAlgebra
using Statistics
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

num_t_passed = 360*4
experiments = [360 for n=1:4]

"""
experiments = [240, 240, 240]
"""

experiment_results_alpha = []
experiment_results_beta = []
experiment_results_gamma = []
experiment_results_objective = []

# iterations out: 4, 6, 9, 10, 11
# infeasible: 6
"""
num_t_passed = 4*240
experiments = [240]
"""

alpha_global = zeros(num_z*num_tech)
beta_global = zeros(num_z*num_tech)
gamma_global = zeros(num_z*num_tech)

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

    g_obs = vcat(albe_obs, alde_obs, at_obs, be_obs, cz_obs, de_obs, fr_obs, hr_obs, hu_obs, nl_obs, pl_obs, ro_obs, si_obs, sk_obs)
    g_obs = convert(Vector{Float64}, g_obs)

    demand = vec(demand_g[(num_t_passed+1):(num_t_passed)+num_t, :])
    lambda_obs = vec(lambda_obs_g[(num_t_passed+1):(num_t_passed)+num_t, :])
    np_obs = vec(np_obs_g[(num_t_passed+1):(num_t_passed)+num_t, :])
    demand = convert(Vector{Float64}, demand)
    lambda_obs = convert(Vector{Float64}, lambda_obs)
    np_obs = convert(Vector{Float64}, np_obs)

    ptdf_z = vec(ptdf_z_g[(num_t_passed*num_j+1):(num_t_passed*num_j)+num_t*num_j, :])
    ram = vec(ram_g[(num_t_passed*num_j+1):(num_t_passed*num_j)+num_t*num_j, :])
    ptdf_z = convert(Vector{Float64}, ptdf_z)
    ram = convert(Vector{Float64}, ram)

    ren_gen = vec(ren_gen_g[(num_t_passed+1):(num_t_passed)+num_t, :])
    ren_gen = convert(Vector{Float64}, ren_gen)

    g_max_t = zeros(num_z*num_tech*num_t)
    for z in 1:num_z
        for t in 1:num_t
            for tech in 1:num_tech
                g_max_t[num_t*num_tech*(z-1)+num_t*(tech-1)+t] = g_max_g[num_tech*(z-1)+tech] 
            end
        end
    end

    #model = direct_model(Gurobi.Optimizer())
    #set_time_limit_sec(model, 120.0)
    #set_optimizer_attribute(model, "NonConvex", 2)
    model = Model(HiGHS.Optimizer)
    set_optimizer_attribute(model, "presolve", "on")
    set_optimizer_attribute(model, "time_limit", 180.0)

    @variable(model, c[1:num_z*num_tech*num_t] >= 0)

    @variable(model, alpha[1:num_z*num_tech])
    @variable(model, beta[1:num_z*num_tech])
    @variable(model, gamma[1:num_z*num_tech])

    @constraint(model, alpha .+ alpha_global .>= 0)
    @constraint(model, beta .+ beta_global .>= 0)
    @constraint(model, gamma .+ gamma_global .>= 0)

    for z in 1:num_z
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

    @variable(model, eps1[1:num_z*num_t] >= 0)
    @variable(model, eps2[1:num_z*num_t] >= 0)
    #@variable(model, eps3[1:num_z*num_t] >= 0)
    #@variable(model, eps4[1:num_z*num_t] >= 0)

    @variable(model, g[1:num_z*num_tech*num_t] >= 0)
    @variable(model, np[1:num_z*num_t])

    # dual variables
    @variable(model, lambda[1:num_z*num_t])
    @variable(model, lambda_exchange[1:num_t])

    @variable(model, mu_gen[1:num_z*num_tech*num_t] <= 0)
    @variable(model, mu_exchange[1:num_j*num_t] <= 0)

    # for minimasing duality gap
    @variable(model, epsilon_duality_abs)
    @variable(model, epsilon_duality)

    @constraint(model, g .== g_obs) # observe generation

    A_balance = spzeros(num_z*num_t, num_z*num_tech*num_t+num_z*num_t) # contains g and np
    prev_pos = num_z*num_tech*num_t
    for z in 1:num_z
        for t in 1:num_t
            for tech in 1:num_tech
                A_balance[num_t*(z-1)+t, num_t*num_tech*(z-1)+num_t*(tech-1)+t] = 1
            end
            A_balance[num_t*(z-1)+t, prev_pos+num_t*(z-1)+t] = -1 # np
        end
    end

    b1_balance = demand - ren_gen

    B_gen = sparse(cat(Matrix(I, num_z*num_tech*num_t, num_z*num_tech*num_t), spzeros(num_z*num_tech*num_t, num_z*num_t); dims=(2)))
    b2_gen = g_max_t

    A_exchange = spzeros(num_t, num_z*num_tech*num_t+num_z*num_t)
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
    B_exchange = sparse(cat(spzeros(num_j*num_t, num_z*num_tech*num_t), B_exchange_temp; dims=(2)))

    # primal constraints
    #@constraint(model, balance, A_balance * vcat(g, np) .== b1_balance) # demand equality
    @constraint(model, B_exchange * vcat(g, np) .<= b2_exchange) # net position, rams and ptdfs
    @constraint(model, sum_z_np(np, num_t) .== 0)

    # dual constraints
    @constraint(model, cat(A_balance, A_exchange; dims=(1))' * vcat(lambda, lambda_exchange) .+ cat(B_gen, B_exchange; dims=(1))' * vcat(mu_gen, mu_exchange) .== vcat(c, spzeros(num_z*num_t)))

    # strong duality gap theorem with epsilon
    @constraint(model, sum(c) .- b1_balance' * lambda .- b1_exchange' * lambda_exchange .- b2_gen' * mu_gen .- b2_exchange' * mu_exchange == epsilon_duality_abs)
    @constraint(model, epsilon_duality_abs <= epsilon_duality)
    @constraint(model, -1*epsilon_duality_abs <= epsilon_duality)

    @constraint(model, lambda .- lambda_obs .== eps1 .- eps2)
    #@constraint(model, np .- np_obs .== eps3 .- eps4)

    u = ones(num_z*num_t)
    @objective(model, Min, eps1' * u + eps2' * u + epsilon_duality)
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
        push!(experiment_results_objective, JuMP.objective_value(model))
    end
   
    global iteration_count += 1  
end

# SAVE COEFFICIENTS

objective_weights = 1 .- normalize(experiment_results_objective)

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

        alpha_mean = mean(alpha_exp[findall(!iszero, alpha_exp)], Weights(objective_weights[findall(!iszero, alpha_exp)]))
        beta_mean = mean(beta_exp[findall(!iszero, beta_exp)], Weights(objective_weights[findall(!iszero, beta_exp)]))
        gamma_mean = mean(gamma_exp[findall(!iszero, gamma_exp)], Weights(objective_weights[findall(!iszero, gamma_exp)]))

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

XLSX.writetable("cost_coefficients_large_winter.xlsx",
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