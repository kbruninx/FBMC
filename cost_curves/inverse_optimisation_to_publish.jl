using JuMP, HiGHS
using DataFrames, XLSX
using LinearAlgebra
using Statistics
using SparseArrays
using StatsBase

# setting up the solver engine
model = Model(HiGHS.Optimizer)
set_optimizer_attribute(model, "presolve", "on")
set_optimizer_attribute(model, "time_limit", 180.0)

# cost curve variables (coefficients)
@variable(model, c[1:(num_z+num_z_non_fbmc)*num_tech*num_t] >= 0)
@variable(model, alpha[1:(num_z+num_z_non_fbmc)*num_tech] >= 0)
@variable(model, beta[1:(num_z+num_z_non_fbmc)*num_tech] >= 0)
@variable(model, gamma[1:(num_z+num_z_non_fbmc)*num_tech] >= 0)

# quadratic curve constraints for different technologies (with different fuel)
for z in 1:(num_z+num_z_non_fbmc)
    for tech in 3:num_tech
        if tech == 7
            @constraint(model, gamma[num_tech*(z-1)+tech] + gamma[num_tech*(z-1)+tech] == 0)
        end
    end
    for t in 1:num_t
        for tech in 1:num_tech
            if tech == 1 || tech == 8 || tech == 9 || tech == 10 # without fuel-responsiveness
                @constraint(model, c[num_t*num_tech*(z-1)+num_t*(tech-1)+t] == alpha[num_tech*(z-1)+tech] + (beta[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t] + gamma[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t]^2))
            elseif tech == 2 || tech == 3 || tech == 5 # coal
                @constraint(model, c[num_t*num_tech*(z-1)+num_t*(tech-1)+t] == coal_prices[t]/10 * (alpha[num_tech*(z-1)+tech] + beta[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t] + gamma[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t]^2))
            elseif tech == 4 # gas
                @constraint(model, c[num_t*num_tech*(z-1)+num_t*(tech-1)+t] == gas_prices[t]/10 * (alpha[num_tech*(z-1)+tech] + beta[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t] + gamma[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t]^2))
            elseif tech == 6 # oil
                @constraint(model, c[num_t*num_tech*(z-1)+num_t*(tech-1)+t] == oil_prices[t]/10 * (alpha[num_tech*(z-1)+tech] + beta[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t] + gamma[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t]^2))
            elseif tech == 7 # hydro
                @constraint(model, c[num_t*num_tech*(z-1)+num_t*(tech-1)+t] == alpha[num_tech*(z-1)+tech] + beta[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t])
            end
        end
    end
end

# for L norms, replacing tha absolute expression in the objective
@variable(model, eps1[1:(num_z+num_z_non_fbmc)*num_t] >= 0)
@variable(model, eps2[1:(num_z+num_z_non_fbmc)*num_t] >= 0)

@variable(model, g[1:(num_z+num_z_non_fbmc)*num_tech*num_t] >= 0)
@variable(model, np[1:num_z*num_t])
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

@constraint(model, g .== g_obs) # observe generation levels

# CONSTRUCTING MATRICES FOR THE INVERSE FORMULATION

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

B_exchange_temp = spzeros(num_j*num_t, num_z*num_t)
ram = spzeros(num_j*num_t)
for t in 1:num_t
    df_ptdf_h = df_ptdf[df_ptdf.DateTime .== timestamps[num_t_passed+t], :]
    for j in 1:size(df_ptdf_h)[1]
        for z in 1:num_z
            B_exchange_temp[num_t*(j-1)+t, num_t*(z-1)+t] = df_ptdf_h[j, fbmc_zones[z]]
            ram[num_t*(j-1)+t] = df_ptdf_h[j, :ram]
        end
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
@constraint(model, balance, A_balance * vcat(g, np) .== b1_balance) # demand equality
@constraint(model, B_exchange * vcat(g, np, atc_ex_1, atc_ex_2) .<= b2_exchange) # net position, rams and ptdfs
@constraint(model, sum_z_np(np, num_t) .== 0)

# dual constraints
@constraint(model, cat(A_balance, A_exchange; dims=(1))' * vcat(lambda, lambda_exchange) .+ cat(B_gen, B_exchange; dims=(1))' * vcat(mu_gen, mu_exchange) .== vcat(c, spzeros(num_z*num_t), spzeros(2*num_atc_border*num_t)))

# relaxed duality
@constraint(model, sum(c) .- b1_balance' * lambda .- b1_exchange' * lambda_exchange .- b2_gen' * mu_gen .- b2_exchange' * mu_exchange == epsilon_duality_abs)
#@constraint(model, epsilon_duality_abs == 0) # strong duality theorem, for L norm only

@constraint(model, epsilon_duality_abs <= epsilon_duality)
@constraint(model, -1*epsilon_duality_abs <= epsilon_duality)
@constraint(model, lambda .- lambda_obs .== eps1 .- eps2) # for L1 norm

u = ones((num_z+num_z_non_fbmc)*num_t)

#@objective(model, Min, sum((lambda - lambda_obs) .^ 2)) # L2 norm only
#@objective(model, Min, eps1' * u + eps2' * u) # L1 norm only
@objective(model, Min, eps1' * u + eps2' * u + epsilon_duality) # L1 norm and duality gap minimisation
#@objective(model, Min, sum((lambda - lambda_obs) .^ 2) + epsilon_duality) # L2 norm and duality gap minimisation

optimize!(model)