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

num_t_passed = 0
num_t = 1

coal_prices = coal_prices_g[(num_t_passed+1):(num_t_passed)+num_t]
oil_prices = oil_prices_g[(num_t_passed+1):(num_t_passed)+num_t]
gas_prices = gas_prices_g[(num_t_passed+1):(num_t_passed)+num_t]
eua_prices = eua_prices_g[(num_t_passed+1):(num_t_passed)+num_t]

demand = vec(demand_g[(num_t_passed+1):(num_t_passed)+num_t, :])
demand = convert(Vector{Float64}, demand)

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

nl_solver = optimizer_with_attributes(Ipopt.Optimizer, "print_level"=>0)
mip_solver = optimizer_with_attributes(Gurobi.Optimizer, "output_flag"=>false)
minlp_solver = optimizer_with_attributes(Juniper.Optimizer, "nl_solver"=>nl_solver, "mip_solver"=>mip_solver)

# MIP optimizer
gurobi = optimizer_with_attributes(Gurobi.Optimizer, 
                                         MOI.Silent() => true,
                                         "Presolve"   => 1) 

# NLP optimizer
ipopt = optimizer_with_attributes(Ipopt.Optimizer, 
                                        MOI.Silent() => true, 
                                        "sb" => "yes", 
                                        "max_iter"   => 9999)

# Global optimizer
alpine = optimizer_with_attributes(Alpine.Optimizer, 
                                         "nlp_solver" => ipopt,
                                         "mip_solver" => gurobi,
                                         "minlp_solver" => minlp_solver)
model = Model(alpine)

@variable(model, c[1:num_z*num_tech*num_t] >= 0)
@variable(model, g[1:num_z*num_tech*num_t] >= 0.0)
@variable(model, u[1:num_z*num_tech*num_t], Bin)

@variable(model, np[1:num_z*num_t])

for z in 1:num_z
    for t in 1:num_t
        for tech in 1:num_tech
            if tech == 2 || tech == 3 || tech == 5 # coal
                @NLconstraint(model, c[num_t*num_tech*(z-1)+num_t*(tech-1)+t] == u[num_t*num_tech*(z-1)+num_t*(tech-1)+t] * (alpha[num_tech*(z-1)+tech] + coal_prices[t]/10 * (beta[num_tech*(z-1)+tech] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t] + gamma[num_tech*(z-1)+tech] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t]^2)))
            elseif tech == 4 # gas
                @NLconstraint(model, c[num_t*num_tech*(z-1)+num_t*(tech-1)+t] == u[num_t*num_tech*(z-1)+num_t*(tech-1)+t] * (alpha[num_tech*(z-1)+tech] + gas_prices[t]/10 * (beta[num_tech*(z-1)+tech] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t] + gamma[num_tech*(z-1)+tech] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t]^2)))
            elseif tech == 6 # oil
                @NLconstraint(model, c[num_t*num_tech*(z-1)+num_t*(tech-1)+t] == u[num_t*num_tech*(z-1)+num_t*(tech-1)+t] * (alpha[num_tech*(z-1)+tech] + oil_prices[t]/10 * (beta[num_tech*(z-1)+tech] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t] + gamma[num_tech*(z-1)+tech] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t]^2)))
            else
                @NLconstraint(model, c[num_t*num_tech*(z-1)+num_t*(tech-1)+t] == u[num_t*num_tech*(z-1)+num_t*(tech-1)+t] * (alpha[num_tech*(z-1)+tech] + (beta[num_tech*(z-1)+tech] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t] + gamma[num_tech*(z-1)+tech] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t]^2)))
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
@constraint(model, cat(B_gen, B_exchange; dims=(1)) * vcat(g, np) .<= vcat(b2_gen .* u, b2_exchange)) # combined inequality constraint

@constraint(model, sum_z_np(np, num_t) .== 0)

@objective(model, Min, sum(c))

optimize!(model)