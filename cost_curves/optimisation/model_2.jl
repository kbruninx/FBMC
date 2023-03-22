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

c_init = zeros(num_z*num_tech*num_t)

"""
merit_order = [7,8,1,9,10,2,5,3,4,6]

for z in 1:num_z
    for t in 1:num_t
        step = 1
        for m in merit_order
            c_init[num_t*num_tech*(z-1)+num_t*(m-1)+t] = 1*step
            step += 1
        end
    end
end
"""

model = Model(Gurobi.Optimizer)
set_optimizer_attribute(model, "NonConvex", 2)

@variable(model, c[1:num_z*num_tech*num_t] >= 0)

# auxiliary variables to represent the norm
@variable(model, alpha[1:num_z*num_tech*num_t] >= 0)
@variable(model, beta[1:num_z*num_tech*num_t] >= 0)
@variable(model, gamma[1:num_z*num_t] >= 0)
@variable(model, delta[1:num_z*num_t] >= 0)

@variable(model, g[1:num_z*num_tech*num_t] >= 0)
@variable(model, np[1:num_z*num_t])

#@constraint(model, g .== g_obs) # observe generation
#@constraint(model, np .== np_obs) # observe net positions

# dual variables
@variable(model, lambda_exchange[1:num_t])
@variable(model, mu_gen[1:num_z*num_tech*num_t] <= 0)
@variable(model, mu_exchange[1:num_j*num_t] <= 0)

"""
function generate_merit_order(t, z, tech_a, tech_b)
    @constraint(model, c[num_t*num_tech*(z-1)+num_t*(tech_a-1)+t] <= c[num_t*num_tech*(z-1)+num_t*(tech_b-1)+t])
    @constraint(model, c[num_t*num_tech*(z-1)+num_t*(tech_b-1)+t] - c[num_t*num_tech*(z-1)+num_t*(tech_a-1)+t] >= 0.1)
end

# merit order
for z in 1:num_z
    for t in 1:num_t
        generate_merit_order(t, z, 8, 1)
        generate_merit_order(t, z, 1, 2)
        generate_merit_order(t, z, 2, 5)
        generate_merit_order(t, z, 5, 3)
        generate_merit_order(t, z, 3, 4)
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
@constraint(model, cat(A_balance, A_exchange; dims=(1))' * vcat(lambda_obs, lambda_exchange) .+ cat(B_gen, B_exchange; dims=(1))' * vcat(mu_gen, mu_exchange) .== vcat(c, spzeros(num_z*num_t)))

# strong duality gap theorem
@constraint(model, c' * g .- b1_balance' * lambda_obs .- b1_exchange' * lambda_exchange .- b2_gen' * mu_gen .- b2_exchange' * mu_exchange == 0)

#@constraint(model, c .- c_init .== alpha .- beta)
@constraint(model, g .- g_obs .== alpha .- beta)
#@constraint(model, np .- np_obs .== gamma .- delta)

# primal constraints
@constraint(model, balance, A_balance * vcat(g, np) .== b1_balance) # demand equality
@constraint(model, B_exchange * vcat(g, np) .<= b2_exchange) # net position, rams and ptdfs
#@constraint(model, cat(B_gen, B_exchange; dims=(1)) * vcat(g, np) .<= vcat(b2_gen, b2_exchange)) # combined inequality constraint
@constraint(model, sum_z_np(np) .== 0)

u = ones(num_z*num_tech*num_t)
#u2 = ones(num_z*num_t)
#@objective(model, Min, alpha' * u + beta' * u + gamma' * u2 + delta' * u2)
@objective(model, Min, alpha' * u + beta' * u )

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