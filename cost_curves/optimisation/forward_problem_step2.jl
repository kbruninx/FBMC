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
using Juniper

model = Model(Ipopt.Optimizer)
#model = Model(Gurobi.Optimizer)
#set_optimizer_attribute(model, "tol", 1e-10)
#set_optimizer_attribute(model, "NonConvex", 2)

@variable(model, c[1:num_z*num_tech*num_t] >= 0)
@variable(model, g[i = 1:num_z*num_tech*num_t] >= 0, start = g_obs[i])
#@variable(model, g[1:num_z*num_tech*num_t] >= 0)

#set_start_value(g, g_prev)

@variable(model, np[1:num_z*num_t])

for z in 1:num_z
    for t in 1:num_t
        for tech in 1:num_tech
            if tech == 2 || tech == 3 || tech == 5 # coal
                @constraint(model, c[num_t*num_tech*(z-1)+num_t*(tech-1)+t] == u[num_t*num_tech*(z-1)+num_t*(tech-1)+t] * (alpha[num_tech*(z-1)+tech] + coal_prices[t]/10 * (beta[num_tech*(z-1)+tech] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t] + gamma[num_tech*(z-1)+tech] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t]^2)))
            elseif tech == 4 # gas
                @constraint(model, c[num_t*num_tech*(z-1)+num_t*(tech-1)+t] == u[num_t*num_tech*(z-1)+num_t*(tech-1)+t] * (alpha[num_tech*(z-1)+tech] + gas_prices[t]/10 * (beta[num_tech*(z-1)+tech] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t] + gamma[num_tech*(z-1)+tech] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t]^2)))
            elseif tech == 6 # oil
                @constraint(model, c[num_t*num_tech*(z-1)+num_t*(tech-1)+t] == u[num_t*num_tech*(z-1)+num_t*(tech-1)+t] * (alpha[num_tech*(z-1)+tech] + oil_prices[t]/10 * (beta[num_tech*(z-1)+tech] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t] + gamma[num_tech*(z-1)+tech] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t]^2)))
            else
                @constraint(model, c[num_t*num_tech*(z-1)+num_t*(tech-1)+t] == u[num_t*num_tech*(z-1)+num_t*(tech-1)+t] * (alpha[num_tech*(z-1)+tech] + (beta[num_tech*(z-1)+tech] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t] + gamma[num_tech*(z-1)+tech] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t]^2)))
            end
        end
    end
end

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
@constraint(model, balance, A_balance * vcat(g, np) .== b1_balance)
#@constraint(model, B_exchange * vcat(g, np) .<= b2_exchange) # net position, rams and ptdfs
@constraint(model, cat(B_gen, B_exchange; dims=(1)) * vcat(g, np) .<= vcat(b2_gen .* u, b2_exchange)) # combined inequality constraint

@constraint(model, sum_z_np(np, num_t) .== 0)

@objective(model, Min, sum(c))
"""
@NLobjective(
    model, 
    Min, 
    sum(
        (
            (
                (
                    tech == 2 || tech == 3 || tech == 5 ? (
                        alpha[num_tech*(z-1)+tech] + coal_prices[t]/10 * (beta[num_tech*(z-1)+tech] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t] + gamma[num_tech*(z-1)+tech] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t]^2)
                    ) : (
                        tech == 4 ? (
                            alpha[num_tech*(z-1)+tech] + gas_prices[t]/10 * (beta[num_tech*(z-1)+tech] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t] + gamma[num_tech*(z-1)+tech] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t]^2)
                        ) : (
                            tech == 6 ? (
                                alpha[num_tech*(z-1)+tech] + oil_prices[t]/10 * (beta[num_tech*(z-1)+tech] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t] + gamma[num_tech*(z-1)+tech] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t]^2)
                            ) : alpha[num_tech*(z-1)+tech] + (beta[num_tech*(z-1)+tech] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t] + gamma[num_tech*(z-1)+tech] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t]^2)
                        )
                    )
                ) for tech in 1:num_tech
            ) for z in 1:num_z
        ) for t in 1:num_t
    )
)
"""

optimize!(model)

np_obs = JuMP.value.(np)
g_obs = JuMP.value.(g)
u = zeros(num_t*num_tech*num_z)
for z in 1:num_z
    for t in 1:num_t
        for tech in 1:num_tech
            if g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t] > 0
                u[num_t*num_tech*(z-1)+num_t*(tech-1)+t] = 1
            end
        end
    end
end