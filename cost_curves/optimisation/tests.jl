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

cost_array = JuMP.value.(c)
costs = []
generation = []
max_gen = []
z = 9
t = 10

for tech in 1:num_tech
    push!(costs, cost_array[num_t*num_tech*(z-1)+num_t*(tech-1)+t])
    push!(generation, JuMP.value.(g)[num_t*num_tech*(z-1)+num_t*(tech-1)+t])
    push!(max_gen, g_max_t[num_t*num_tech*(z-1)+num_t*(tech-1)+t])
end

df_offer = DataFrames.DataFrame(cost=costs, gen=generation, max_gen=max_gen)
sort!(df_offer, [:cost])

demand2 = b1_balance[num_t*(z-1)+t]

println(df_offer)
println(demand2)