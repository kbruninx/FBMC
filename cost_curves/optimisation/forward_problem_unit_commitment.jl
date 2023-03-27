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

model = Model(Gurobi.Optimizer)

@variable(model, c[1:num_z*num_tech*num_t] >= 0)
@variable(model, g_obs[i]-0.1 <= g[i = 1:num_z*num_tech*num_t] <= u[i] * (g_obs[i]+0.1))
@constraint(model, g .>= 0)

@variable(model, np[1:num_z*num_t])
@constraint(model, np .== np_obs)
#@constraint(model, g .== g_obs)
@constraint(model, c .== c_obs .* u)

"""
for z in 1:num_z
    for t in 1:num_t
        for tech in 1:num_tech
            if tech == 2 || tech == 3 || tech == 5 # coal
                @constraint(model, c[num_t*num_tech*(z-1)+num_t*(tech-1)+t] == alpha[num_tech*(z-1)+tech] + coal_prices[t]/10 * (beta[num_tech*(z-1)+tech] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t] + gamma[num_tech*(z-1)+tech] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t]^2))
            elseif tech == 4 # gas
                @constraint(model, c[num_t*num_tech*(z-1)+num_t*(tech-1)+t] == alpha[num_tech*(z-1)+tech] + gas_prices[t]/10 * (beta[num_tech*(z-1)+tech] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t] + gamma[num_tech*(z-1)+tech] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t]^2))
            elseif tech == 6 # oil
                @constraint(model, c[num_t*num_tech*(z-1)+num_t*(tech-1)+t] == alpha[num_tech*(z-1)+tech] + oil_prices[t]/10 * (beta[num_tech*(z-1)+tech] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t] + gamma[num_tech*(z-1)+tech] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t]^2))
            else
                @constraint(model, c[num_t*num_tech*(z-1)+num_t*(tech-1)+t] == alpha[num_tech*(z-1)+tech] + (beta[num_tech*(z-1)+tech] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t] + gamma[num_tech*(z-1)+tech] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t]^2))
            end
        end
    end
end
"""

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

B_gen = sparse(Matrix(I, num_z*num_tech*num_t, num_z*num_tech*num_t))
b2_gen = g_max_t

# primal constraints
@constraint(model, balance, A_balance * vcat(g, np) .== b1_balance)
#@constraint(model, B_gen * g .<= u .* b2_gen) # combined inequality constraint

#@objective(model, Min, sum(c))

"""
@objective(
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

"""
@objective(
    model, 
    Min, 
    sum(
        (
            tech == 2 || tech == 3 || tech == 5 ? (
                u[num_t*num_tech*(z-1)+num_t*(tech-1)+t] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t] * coal_prices[t]/10 * (beta[num_tech*(z-1)+tech] + 2 * gamma[num_tech*(z-1)+tech] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t])
            ) : (
                tech == 4 ? (
                    u[num_t*num_tech*(z-1)+num_t*(tech-1)+t] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t] * gas_prices[t]/10 * (beta[num_tech*(z-1)+tech] + 2 * gamma[num_tech*(z-1)+tech] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t])
                ) : (
                    tech == 6 ? (
                        u[num_t*num_tech*(z-1)+num_t*(tech-1)+t] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t] * oil_prices[t]/10 * (beta[num_tech*(z-1)+tech] + 2 * gamma[num_tech*(z-1)+tech] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t])
                    ) : u[num_t*num_tech*(z-1)+num_t*(tech-1)+t] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t] * (beta[num_tech*(z-1)+tech] + 2 * gamma[num_tech*(z-1)+tech] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t])
                )
            )) for tech in 1:num_tech for z in 1:num_z for t in 1:num_t
    )
)
"""

@objective(
    model, 
    Min, 
    c' * g
)

optimize!(model)