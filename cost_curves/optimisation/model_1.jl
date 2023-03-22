using JuMP, BilevelJuMP, Gurobi, Dualization
using XLSX
using LinearAlgebra
using Alpine
using Ipopt
using Statistics
using QuadraticToBinary
using Plots
using SparseArrays
using Formatting

c_init = zeros(num_z*num_tech*num_t)

model = Model(Gurobi.Optimizer)

@variable(model, c[1:num_z*num_tech*num_t] >= 0)
@variable(model, g[1:num_z*num_tech*num_t] >= 0)
@variable(model, np[1:num_z*num_t])

@variable(model, eps1[1:num_z*num_tech*num_t] >= 0)
@variable(model, eps2[1:num_z*num_tech*num_t] >= 0)

@variable(model, alpha[1:num_z*num_tech] >= 0)
@variable(model, beta[1:num_z*num_tech] >= 0)
@variable(model, gamma[1:num_z*num_tech] >= 0)

# dual variables
@variable(model, lambda_exchange[1:num_t])
@variable(model, mu_gen[1:num_z*num_tech*num_t] <= 0)
@variable(model, mu_exchange[1:num_j*num_t] <= 0)

# for minimasing duality gap
@variable(model, epsilon_abs)
@variable(model, epsilon)

@constraint(model, g .== g_obs)
#@constraint(model, np .== np_obs)

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
@constraint(model, balance, A_balance * vcat(g_obs, np) .== b1_balance)

B_gen = sparse(cat(Matrix(I, num_z*num_tech*num_t, num_z*num_tech*num_t), spzeros(num_z*num_tech*num_t, num_z*num_t); dims=(2)))
b2_gen = g_max_t
#@constraint(model, B_gen * vcat(g, np) .<= b2_gen) # maximum generation capacity for conventional plants

A_exchange = spzeros(num_t, num_z*num_tech*num_t+num_z*num_t)
for t in 1:num_t
    for z in 1:num_z
        A_exchange[t, prev_pos+num_t*(z-1)+t] = 1
    end
end

b1_exchange = spzeros(num_t)
#@constraint(model, A_exchange * np .== b1_exchange)

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

@constraint(model, B_exchange * vcat(g, np) .<= b2_exchange) # net position, rams and ptdfs

for z in 1:num_z
    for tech in 3:num_tech
        if tech == 7
            @constraint(model, beta[num_tech*(z-1)+tech] == 0)
            @constraint(model, gamma[num_tech*(z-1)+tech] == 0)
        end
    end
    for t in 1:num_t
        for tech in 1:num_tech
            if tech == 1 || tech == 8 || tech == 9 || tech == 10 # without fuel-responsiveness
                @constraint(model, c[num_t*num_tech*(z-1)+num_t*(tech-1)+t] == alpha[num_tech*(z-1)+tech] + (beta[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t] + gamma[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t]^2))
            elseif tech == 2 || tech == 3 || tech == 5 # coal
                @constraint(model, c[num_t*num_tech*(z-1)+num_t*(tech-1)+t] == alpha[num_tech*(z-1)+tech] + coal_prices[t]/100 * (beta[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t] + gamma[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t]^2))
            elseif tech == 4 # gas
                @constraint(model, c[num_t*num_tech*(z-1)+num_t*(tech-1)+t] == alpha[num_tech*(z-1)+tech] + gas_prices[t]/100 * (beta[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t] + gamma[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t]^2))
            elseif tech == 6 # oil
                @constraint(model, c[num_t*num_tech*(z-1)+num_t*(tech-1)+t] == alpha[num_tech*(z-1)+tech] + oil_prices[t]/100 * (beta[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t] + gamma[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t]^2))
            else # tech == 7 # hydro
                @constraint(model, c[num_t*num_tech*(z-1)+num_t*(tech-1)+t] == alpha[num_tech*(z-1)+tech])
            end
        end
    end
end

# dual feasibility
@constraint(model, cat(A_balance, A_exchange; dims=(1))' * vcat(lambda_obs, lambda_exchange) .+ cat(B_gen, B_exchange; dims=(1))' * vcat(mu_gen, mu_exchange) .<= vcat(c, spzeros(num_z*num_t))) # for generation
# for np they are already met -> constant <= 0
#@constraint(model, A_exchange' * lambda_exchange .+ B_exchange' * mu_exchange .== 0) # for np


# duality gap
@constraint(model, sum(c) .- b1_balance' * lambda_obs .- b1_exchange' * lambda_exchange .- b2_gen' * mu_gen .- b2_exchange' * mu_exchange == 0)

"""
@constraint(model, epsilon_abs <= epsilon)
@constraint(model, -1*epsilon_abs <= epsilon)
@objective(model, Min, epsilon)
"""

@constraint(model, c .- c_init .== eps1 .- eps2)
u = ones(num_z*num_tech*num_t)
@objective(model, Min, eps1' * u + eps2' * u )


optimize!(model)

function show_stuff(z, t)
    bids = []
    for tech in 1:num_tech
        amount = g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t]
        cost = JuMP.value.(alpha)[num_tech*(z-1)+tech] + fuel_price_oil * (JuMP.value.(beta)[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t] + JuMP.value.(gamma)[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t]^2)
        push!(bids, (amount, cost))
    end

    sort!(bids, by = x -> x[2])

    x = zeros(num_tech)
    y = zeros(num_tech)

    for tech in 1:num_tech
        if tech > 1
            x[tech] = x[tech-1]+bids[tech][1]
        else
            x[tech] = bids[tech][1]
        end
        y[tech] = bids[tech][2]
    end

    plot(x,y)
end

function show_cost_curve(z, tech)
    x = range(0, round(g_max[num_tech*(z-1)+tech]), length=100)
    y = JuMP.value.(alpha)[num_tech*(z-1)+tech] .+ JuMP.value.(beta)[num_tech*(z-1)+tech] .* x .+ JuMP.value.(gamma)[num_tech*(z-1)+tech] .* x .^ 2
    plot(x, y)
end