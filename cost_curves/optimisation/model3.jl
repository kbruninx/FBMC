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

zones =  ["ALBE", "ALDE", "AT", "BE", "CZ", "DE_LU", "FR", "HR", "HU", "NL", "PL", "RO", "SI", "SK"]

plant_types = [
    "biomass", "brown_coal", "coal_gas", "natural_gas", "hard_coal", "oil", "hydro", 
    "nuclear", "waste", "other"
] 

function sum_z_np(np)
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

model = Model(Gurobi.Optimizer)
set_optimizer_attribute(model, "NonConvex", 2)

@variable(model, c[1:num_z*num_tech*num_t] >= 0)

@variable(model, alpha[1:num_z*num_tech] >= 0)
@variable(model, beta[1:num_z*num_tech] >= 0)
@variable(model, gamma[1:num_z*num_tech] >= 0)

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
                @constraint(model, c[num_t*num_tech*(z-1)+num_t*(tech-1)+t] == alpha[num_tech*(z-1)+tech] + coal_prices[t]/10 * (beta[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t] + gamma[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t]^2))
            elseif tech == 4 # gas
                @constraint(model, c[num_t*num_tech*(z-1)+num_t*(tech-1)+t] == alpha[num_tech*(z-1)+tech] + gas_prices[t]/10 * (beta[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t] + gamma[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t]^2))
            elseif tech == 6 # oil
                @constraint(model, c[num_t*num_tech*(z-1)+num_t*(tech-1)+t] == alpha[num_tech*(z-1)+tech] + oil_prices[t]/10 * (beta[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t] + gamma[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t]^2))
            else # tech == 7 # hydro
                @constraint(model, c[num_t*num_tech*(z-1)+num_t*(tech-1)+t] == alpha[num_tech*(z-1)+tech])
            end
        end
    end
end

# auxiliary variables to represent the norm
#@variable(model, alpha[1:num_z*num_tech*num_t] >= 0)
#@variable(model, beta[1:num_z*num_tech*num_t] >= 0)
#@variable(model, gamma[1:num_z*num_tech*num_t] >= 0)
#@variable(model, delta[1:num_z*num_tech*num_t] >= 0)
@variable(model, khi[1:num_z*num_t] >= 0)
@variable(model, ksi[1:num_z*num_t] >= 0)

@variable(model, g[1:num_z*num_tech*num_t] >= 0)
@variable(model, np[1:num_z*num_t])

# dual variables
@variable(model, lambda[1:num_z*num_t])
@variable(model, lambda_exchange[1:num_t])
@variable(model, mu_gen[1:num_z*num_tech*num_t] <= 0)
@variable(model, mu_exchange[1:num_j*num_t] <= 0)

@constraint(model, g .== g_obs) # observe generation
#@constraint(model, lambda .== lambda_obs) # observe prices
#@constraint(model, np .== np_obs) # observe net positions

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

# dual constraints
@constraint(model, cat(A_balance, A_exchange; dims=(1))' * vcat(lambda, lambda_exchange) .+ cat(B_gen, B_exchange; dims=(1))' * vcat(mu_gen, mu_exchange) .== vcat(c, spzeros(num_z*num_t)))

# strong duality gap theorem
@constraint(model, sum(c) .- b1_balance' * lambda .- b1_exchange' * lambda_exchange .- b2_gen' * mu_gen .- b2_exchange' * mu_exchange == 0)

#@constraint(model, c .- c_init .== gamma .- delta)
#@constraint(model, g .- g_obs .== alpha .- beta)
@constraint(model, lambda .- lambda_obs .== khi .- ksi)

# primal constraints
#@constraint(model, balance, A_balance * vcat(g, np) .== b1_balance) # demand equality
@constraint(model, B_exchange * vcat(g, np) .<= b2_exchange) # net position, rams and ptdfs
#@constraint(model, cat(B_gen, B_exchange; dims=(1)) * vcat(g, np) .<= vcat(b2_gen, b2_exchange)) # combined inequality constraint
@constraint(model, sum_z_np(np) .== 0)

u = ones(num_z*num_tech*num_t)
u2 = ones(num_z*num_t)

@objective(model, Min, khi' * u2 + ksi' * u2)
#@objective(model, Min, alpha' * u + beta' * u + gamma' * u + delta' * u)
#@objective(model, Min, alpha' * u + beta' * u )
#@objective(model, Min, gamma' * u + delta' * u )

optimize!(model)


XLSX.openxlsx("marginal_prices_output.xlsx", mode="w") do xf
    for z in 3:num_z
        if z == 3
            sheet = xf[1]
            XLSX.rename!(sheet, zones[z])
        else
            sheet = XLSX.addsheet!(xf, zones[z])  
        end

        sheet["A1"] = "date_time"
        sheet["B1"] = "marginal_type"
        sheet["C1"] = "marginal_cost"
        sheet["D1"] = "production_level"

        i = 1
        for t in 1:num_t
            for tech in 1:num_tech
                if JuMP.value.(c)[num_t*num_tech*(z-1)+num_t*(tech-1)+t] > 10 && JuMP.value.(g)[num_t*num_tech*(z-1)+num_t*(tech-1)+t] > 0
                    sheet[sprintf1("A%d", i+1)] = t
                    sheet[sprintf1("B%d", i+1)] = plant_types[tech]
                    sheet[sprintf1("C%d", i+1)] = JuMP.value.(c)[num_t*num_tech*(z-1)+num_t*(tech-1)+t]
                    sheet[sprintf1("D%d", i+1)] = JuMP.value.(g)[num_t*num_tech*(z-1)+num_t*(tech-1)+t]
                    i += 1
                end
            end
        end
    end
end