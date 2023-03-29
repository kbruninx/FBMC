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

num_t_passed = 0
num_t = 24
day_count = 56

np_obs = zeros(num_z*num_t)
#np_obs = vec(np_obs_g[(num_t_passed+1):(num_t_passed)+num_t, :])
ptdf_z = vec(ptdf_z_g[(num_t_passed*num_j+1):(num_t_passed*num_j)+num_t*num_j, :])
ram = vec(ram_g[(num_t_passed*num_j+1):(num_t_passed*num_j)+num_t*num_j, :])
ptdf_z = convert(Vector{Float64}, ptdf_z)
ram = convert(Vector{Float64}, ram)

# FEASIBILITY PROBLEM

model = Model(Gurobi.Optimizer)
#set_optimizer_attribute(model, "NonConvex", 2)

H = 0.5

@variable(model, np_min[1:num_z*num_t] <= 0)
@variable(model, np_max[1:num_z*num_t])

@variable(model, np_min_eps_plus[1:num_z*num_t] >= 0)
@variable(model, np_min_eps_min[1:num_z*num_t] >= 0)

@variable(model, np_max_eps_plus[1:num_z*num_t] >= 0)
@variable(model, np_max_eps_min[1:num_z*num_t] >= 0)

@constraint(model, np_max .- np_obs .== np_max_eps_plus .- np_max_eps_min)
@constraint(model, np_obs .- np_min .== np_min_eps_plus .- np_min_eps_min)

@constraint(model, np_min .<= np_max)

B_exchange = spzeros(num_j*num_t, num_z*num_t)
b2_exchange = ram
for t in 1:num_t
    for j in 1:num_j
        for z in 1:num_z
            B_exchange[num_t*(j-1)+t, num_t*(z-1)+t] = ptdf_z[num_j*num_t*(z-1) + num_j*(t-1) + j]
        end
    end
end

@constraint(model, B_exchange * np_min .<= b2_exchange)
@constraint(model, B_exchange * np_max .<= b2_exchange)

#@objective(model, Max, sum(H.*(np_max_eps_min .+ np_min_eps_min)) + sum((1-H).*(np_max_eps_plus .+ np_min_eps_plus)))
@objective(model, Max, sum(np_max_eps_plus .+ np_max_eps_min))

optimize!(model)

# OPTIMALITY PROBLEM