using Dates
using JLD
using JuMP, HiGHS, Ipopt, Juniper
using DataFrames, XLSX
using LinearAlgebra
using Statistics
using SparseArrays
using Formatting
using StatsBase

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

function has_interconnector_ex_1(z, b)
    if all_zones[z] == borders[b][1]
        return -1
    elseif all_zones[z] == borders[b][2]
            return 1
    else
        return 0
    end
end

function has_interconnector_ex_2(z, b)
    if all_zones[z] == reverse(borders[b])[1]
        return -1
    elseif all_zones[z] == reverse(borders[b])[2]
            return 1
    else
        return 0
    end
end

function get_atc_ex_1(b, t)
    return df_atc[!, borders[b][1] * "_" * borders[b][2]][t]
end

function get_atc_ex_2(b, t)
    return df_atc[!, borders[b][2] * "_" * borders[b][1]][t]
end

k_closest = 5

num_t_passed = 0
experiments = [240 for n=1:21]

coefficients_data = load(string("coefficients_norm_1_w_atc.jld"))["data"]
experiment_results_alpha = coefficients_data["alpha"]
experiment_results_beta = coefficients_data["beta"]
experiment_results_gamma = coefficients_data["gamma"]
experiment_results_objective = coefficients_data["objective"]
experiment_results_timestamp = coefficients_data["timestamps"]
objective_weights = 1 .- normalize(experiment_results_objective)

alpha = zeros(2*num_tech)
beta = zeros(2*num_tech)
gamma = zeros(2*num_tech)
for z in 3:(num_z+num_z_non_fbmc)
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
    
    alpha = vcat(alpha, alpha_coeffs)
    beta = vcat(beta, beta_coeffs)
    gamma = vcat(gamma, gamma_coeffs)
end

experiment_results_ram_cap = []
experiment_results_fref_demand = []
experiment_results_fref_ren = []
experiment_results_fref_prodout = []
experiment_results_fref_transout = []
experiment_results_objective = []
experiment_results_timestamp = []

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

    g_obs = vcat(albe_obs, alde_obs, at_obs, be_obs, cz_obs, de_obs, fr_obs, hr_obs, hu_obs, nl_obs, pl_obs, ro_obs, si_obs, sk_obs, ch_obs, gb_obs, es_obs, it_obs)
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

    ren_gen_fbmc = vec(ren_gen_g[(num_t_passed+1):(num_t_passed)+num_t, :])
    ren_gen_non_fbmc = vec(ren_gen_g_non_fbmc[(num_t_passed+1):(num_t_passed)+num_t, :])
    ren_gen = vcat(ren_gen_fbmc, ren_gen_non_fbmc)

    ren_gen = convert(Vector{Float64}, ren_gen)

    model = Model(HiGHS.Optimizer)
    set_optimizer_attribute(model, "presolve", "on")
    set_optimizer_attribute(model, "time_limit", 180.0)

    c = Array{Float64}(undef, (num_z+num_z_non_fbmc)*num_tech*num_t) 

    for z in 1:(num_z+num_z_non_fbmc)
        for t in 1:num_t
            for tech in 1:num_tech
                if tech == 1 || tech == 8 || tech == 9 || tech == 10 # without fuel-responsiveness
                    c[num_t*num_tech*(z-1)+num_t*(tech-1)+t] = alpha[num_tech*(z-1)+tech] + (beta[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t] + gamma[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t]^2)
                elseif tech == 2 || tech == 3 || tech == 5 # coal
                    c[num_t*num_tech*(z-1)+num_t*(tech-1)+t] = coal_prices[t]/10 * (alpha[num_tech*(z-1)+tech] + beta[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t] + gamma[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t]^2)
                elseif tech == 4 # gas
                   c[num_t*num_tech*(z-1)+num_t*(tech-1)+t] = gas_prices[t]/10 * (alpha[num_tech*(z-1)+tech] + beta[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t] + gamma[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t]^2)
                elseif tech == 6 # oil
                    c[num_t*num_tech*(z-1)+num_t*(tech-1)+t] = oil_prices[t]/10 * (alpha[num_tech*(z-1)+tech] + beta[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t] + gamma[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t]^2)
                elseif tech == 7 # hydro
                    c[num_t*num_tech*(z-1)+num_t*(tech-1)+t] = alpha[num_tech*(z-1)+tech] + beta[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t]
                end
            end
        end
    end

    @variable(model, eps1_lambda[1:(num_z+num_z_non_fbmc)*num_t] >= 0)
    @variable(model, eps2_lambda[1:(num_z+num_z_non_fbmc)*num_t] >= 0)
    @variable(model, eps1_np[1:num_z*num_t] >= 0)
    @variable(model, eps2_np[1:num_z*num_t] >= 0)
    @variable(model, eps1_ram[1:num_j*num_t] >= 0)
    @variable(model, eps2_ram[1:num_j*num_t] >= 0)

    @variable(model, g[1:(num_z+num_z_non_fbmc)*num_tech*num_t] >= 0)
    @variable(model, np[1:num_z*num_t])

    @variable(model, ram[1:num_j*num_t])
    @variable(model, ram_cap_coeff[1:num_l])
    @variable(model, fref_demand_coeff[1:num_l, 1:num_z])
    @variable(model, fref_ren_coeff[1:num_l, 1:num_z])
    @variable(model, fref_prodout_coeff[1:num_l, 1:k_closest])
    @variable(model, fref_transout_coeff[1:num_l, 1:k_closest])

    @variable(model, atc_ex_1[1:num_atc_border*num_t] >= 0)
    @variable(model, atc_ex_2[1:num_atc_border*num_t] >= 0)

    # dual variables
    @variable(model, lambda[1:(num_z+num_z_non_fbmc)*num_t])
    @variable(model, lambda_exchange[1:num_t])

    @variable(model, mu_gen[1:(num_z+num_z_non_fbmc)*num_tech*num_t] <= 0)
    @variable(model, mu_exchange[1:(num_j*num_t+2*num_atc_border*num_t)] <= 0)

    # for minimasing duality gap
    @variable(model, epsilon_duality_abs)
    @variable(model, epsilon_duality)

    @constraint(model, g .== g_obs) # observe generation

    A_balance = spzeros((num_z+num_z_non_fbmc)*num_t, (num_z+num_z_non_fbmc)*num_tech*num_t+num_z*num_t+2*num_atc_border*num_t) # contains g and np and atc (e/i)
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

    for z in 1:(num_z+num_z_non_fbmc)
        for b in 1:num_atc_border
            for t in 1:num_t
                A_balance[num_t*(z-1)+t, prev_pos+num_t*(b-1)+t] = has_interconnector_ex_1(z, b) # atc exchange direction 1
            end
        end
    end

    prev_pos += num_atc_border*num_t

    for z in 1:(num_z+num_z_non_fbmc)
        for b in 1:num_atc_border
            for t in 1:num_t
                A_balance[num_t*(z-1)+t, prev_pos+num_t*(b-1)+t] = has_interconnector_ex_2(z, b) # atc exchange direction 2
            end
        end
    end

    b1_balance = demand - ren_gen

    B_gen = sparse(cat(Matrix(I, (num_z+num_z_non_fbmc)*num_tech*num_t, (num_z+num_z_non_fbmc)*num_tech*num_t), spzeros((num_z+num_z_non_fbmc)*num_tech*num_t, num_z*num_t + 2*num_atc_border*num_t); dims=(2)))
    b2_gen = g_max_t

    A_exchange = spzeros(num_t, (num_z+num_z_non_fbmc)*num_tech*num_t+num_z*num_t+2*num_atc_border*num_t)

    prev_pos = (num_z+num_z_non_fbmc)*num_tech*num_t
    for t in 1:num_t
        for z in 1:num_z
            A_exchange[t, prev_pos+num_t*(z-1)+t] = 1
        end
    end

    b1_exchange = spzeros(num_t)
    ram_obs = zeros(num_j*num_t)

    B_exchange_temp = spzeros(num_j*num_t, num_z*num_t)
    ram = spzeros(num_j*num_t)
    for t in 1:num_t
        df_ptdf_h = df_ptdf[df_ptdf.DateTime .== timestamps[num_t_passed+t], :]
        for j in 1:size(df_ptdf_h)[1]

            demand_t = zeros(num_z)
            ren_gen_t = zeros(num_z)

            for z in 1:num_z
                B_exchange_temp[num_t*(j-1)+t, num_t*(z-1)+t] = df_ptdf_h[j, fbmc_zones[z]]

                demand_t[z] = demand[num_t*(z-1)+t]
                ren_gen_t[z] = ren_gen[num_t*(z-1)+t]
            end

            l_i = findfirst(item -> item == df_ptdf_h[j, :line_id], observed_lines)
            fref_j_t = fref_demand_coeff[l_i, :]' * demand_t + fref_ren_coeff[l_i, :]' * ren_gen_t + fref_prodout_coeff[l_i, :]' * df_ptdf_k_h[j, ["p"*string(p) for p=1:k_closest]] + fref_transout_coeff[l_i, :]' * df_ptdf_k_h[j, ["t"*string(t) for t=1:(k_closest-1)]]
            ram[num_t*(j-1)+t] = ram_cap_coeff[l_i] * df_ptdf_h[j, :fmax_calc] - df_ptdf_h[j, :frm] - fref_j_t
            ram_obs[num_t*(j-1)+t] = df_ptdf_h[j, :ram]
        end
    end
    
    ram = convert(Vector{Float64}, ram)
    B_exchange_ptdf = sparse(cat(spzeros(num_j*num_t, (num_z+num_z_non_fbmc)*num_tech*num_t), B_exchange_temp, spzeros(num_j*num_t, 2*num_atc_border*num_t); dims=(2)))
    
    B_exchange_atc = spzeros(2*num_atc_border*num_t, (num_z+num_z_non_fbmc)*num_tech*num_t + num_z*num_t + 2*num_atc_border*num_t)
    prev_pos = (num_z+num_z_non_fbmc)*num_tech*num_t + num_z*num_t
    
    atc_1 = spzeros(num_atc_border*num_t)
    for b in 1:num_atc_border
        for t in 1:num_t
            B_exchange_atc[num_t*(b-1)+t, prev_pos+num_t*(b-1)+t] = 1 # atc exchange direction 1
            atc_1[num_t*(b-1)+t] = get_atc_ex_1(b, t)
        end
    end

    prev_pos += num_atc_border*num_t

    atc_2 = spzeros(num_atc_border*num_t)
    for b in 1:num_atc_border
        for t in 1:num_t
            B_exchange_atc[num_t*(b-1)+t, prev_pos+num_t*(b-1)+t] = 1 # atc exchange direction 2
            atc_2[num_t*(b-1)+t] = get_atc_ex_2(b, t)
        end
    end

    b2_exchange = vcat(ram, remove_missing(atc_1), remove_missing(atc_2))
    B_exchange = sparse(cat(B_exchange_ptdf, B_exchange_atc; dims=(1)))

    # primal constraints
    #@constraint(model, balance, A_balance * vcat(g, np) .== b1_balance) # demand equality
    @constraint(model, B_exchange * vcat(g, np, atc_ex_1, atc_ex_2) .<= b2_exchange) # net position, rams and ptdfs
    @constraint(model, sum_z_np(np, num_t) .== 0)

    # dual constraints
    @constraint(model, cat(A_balance, A_exchange; dims=(1))' * vcat(lambda, lambda_exchange) .+ cat(B_gen, B_exchange; dims=(1))' * vcat(mu_gen, mu_exchange) .== vcat(c, spzeros(num_z*num_t), spzeros(2*num_atc_border*num_t)))

    # strong duality gap theorem
    @constraint(model, epsilon_duality_abs == 0) # for L norm only
    @constraint(model, sum(c) .- b1_balance' * lambda .- b1_exchange' * lambda_exchange .- b2_gen' * mu_gen .- b2_exchange' * mu_exchange == epsilon_duality_abs)
    @constraint(model, epsilon_duality_abs <= epsilon_duality)
    @constraint(model, -1*epsilon_duality_abs <= epsilon_duality)

    @constraint(model, lambda .- lambda_obs .== eps1_lambda .- eps2_lambda) # for L1 norm (lambda)
    @constraint(model, np .- np_obs .== eps1_np .- eps2_np) # for L1 norm (np)
    @constraint(model, ram .- ram_obs .== eps1_ram .- eps2_ram) # for L1 norm (ram)

    u_lambda = ones((num_z+num_z_non_fbmc)*num_t)
    u_np = ones(num_z*num_t)
    u_ram = ones(num_j*num_t)

    h_lambda = 1
    h_np = 1
    h_ram = 1

    lambda_obj = h_lambda * (eps1_lambda' * u_lambda + eps2_lambda' * u_lambda)
    np_obj = h_np * (eps1_np' * u_np + eps2_np' * u_np)
    ram_obj = h_ram * (eps1_ram' * u_ram + eps2_ram' * u_ram)
 
    @objective(model, Min, lambda_obj + np_obj + ram_obj) # L1 norm only
    #@objective(model, Min, lambda_obj + np_obj + ram_obj + epsilon_duality) # L1 norm and duality gap minimisation

    # iteratively adjust upon the previous values
    global num_t_passed += exp_len 

    try
        optimize!(model)

        if termination_status(model) == OPTIMAL || termination_status(model) == LOCALLY_SOLVED || termination_status(model) == ALMOST_OPTIMAL || termination_status(model) == ALMOST_LOCALLY_SOLVED
            println("ITERATION ", iteration_count, " FOUND OPTIMAL SOLUTION")
            push!(experiment_results_ram_cap, JuMP.value.(ram_cap_coeff))
            push!(experiment_results_fref_demand, JuMP.value.(fref_demand_coeff))
            push!(experiment_results_fref_ren, JuMP.value.(fref_ren_coeff))
            push!(experiment_results_fref_prodout, JuMP.value.(fref_prodout_coeff))
            push!(experiment_results_fref_transout, JuMP.value.(fref_transout_coeff))
            push!(experiment_results_objective, JuMP.objective_value(model))
            push!(experiment_results_timestamp, timestamps[num_t_passed])
        end
    catch e
        println("ERROR ENCOUNTERED")
    end
   
    global iteration_count += 1  
end

ram_coefficients_data = Dict(
    "observed_lines" => observed_lines,
    "ram_cap_coeff" => experiment_results_ram_cap,
    "fref_demand_coeff" => experiment_results_fref_demand,
    "fref_ren_coeff" => experiment_results_fref_ren,
    "fref_prodout_coeff" => experiment_results_fref_prodout,
    "fref_transout_coeff" => experiment_results_fref_transout,
    "objective" => experiment_results_objective,
    "timestamps" => experiment_results_timestamp
)

save("ram_coefficients.jld", "data", ram_coefficients_data)
