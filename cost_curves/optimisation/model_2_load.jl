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

function remove_missing(mx)
    mx[ismissing.(mx)] .= 0.0
    return mx
end

xf_atc = XLSX.readxlsx("./data/training/atc.xlsx")
xf_da_prices = XLSX.readxlsx("./data/training/day_ahead_prices.xlsx")
xf_da_prices_non_fbmc = XLSX.readxlsx("./data/training/day_ahead_prices_non_fbm.xlsx")
xf_demand = XLSX.readxlsx("./data/training/demand.xlsx")
xf_fuel_prices = XLSX.readxlsx("./data/training/fuel_prices.xlsx")
xf_generation = XLSX.readxlsx("./data/training/generation.xlsx")
xf_ren_gen = XLSX.readxlsx("./data/training/renewable_generation.xlsx")
xf_netpos = XLSX.readxlsx("./data/training/net_positions.xlsx")
xf_ptdf = XLSX.readxlsx("./data/training/ptdfs.xlsx")
xf_capacities = XLSX.readxlsx("./data/installed_capacities_corrected.xlsx")

num_z = 14 # including ALBE and ALDE
num_t = 24 #2906 # size of training set
num_tech = 10
num_j = 133 # maximum amont on CNEs at a given time in the training dataset

coal_prices = vec(remove_missing(xf_fuel_prices["Sheet1"][sprintf1("B2:B%d", num_t+1)]))
oil_prices = vec(remove_missing(xf_fuel_prices["Sheet1"][sprintf1("C2:C%d", num_t+1)]))
gas_prices = vec(remove_missing(xf_fuel_prices["Sheet1"][sprintf1("D2:D%d", num_t+1)]))
eua_prices = vec(remove_missing(xf_fuel_prices["Sheet1"][sprintf1("E2:E%d", num_t+1)]))

coal_prices = convert(Vector{Float64}, coal_prices)
oil_prices = convert(Vector{Float64}, oil_prices)
gas_prices = convert(Vector{Float64}, gas_prices)
eua_prices = convert(Vector{Float64}, eua_prices)

albe_obs = zeros(num_t*num_tech)
alde_obs = zeros(num_t*num_tech)
at_obs = vec(remove_missing(xf_generation["AT"][sprintf1("B2:K%d", num_t+1)]))
be_obs = vec(remove_missing(xf_generation["BE"][sprintf1("B2:K%d", num_t+1)]))
cz_obs = vec(remove_missing(xf_generation["CZ"][sprintf1("B2:K%d", num_t+1)]))
de_obs = vec(remove_missing(xf_generation["DE_LU"][sprintf1("B2:K%d", num_t+1)]))
fr_obs = vec(remove_missing(xf_generation["FR"][sprintf1("B2:K%d", num_t+1)]))
hr_obs = vec(remove_missing(xf_generation["HR"][sprintf1("B2:K%d", num_t+1)]))
hu_obs = vec(remove_missing(xf_generation["HU"][sprintf1("B2:K%d", num_t+1)]))
nl_obs = vec(remove_missing(xf_generation["NL"][sprintf1("B2:K%d", num_t+1)]))
pl_obs = vec(remove_missing(xf_generation["PL"][sprintf1("B2:K%d", num_t+1)]))
ro_obs = vec(remove_missing(xf_generation["RO"][sprintf1("B2:K%d", num_t+1)]))
si_obs = vec(remove_missing(xf_generation["SI"][sprintf1("B2:K%d", num_t+1)]))
sk_obs = vec(remove_missing(xf_generation["SK"][sprintf1("B2:K%d", num_t+1)]))
g_obs = vcat(albe_obs, alde_obs, at_obs, be_obs, cz_obs, de_obs, fr_obs, hr_obs, hu_obs, nl_obs, pl_obs, ro_obs, si_obs, sk_obs)
g_obs = convert(Vector{Float64}, g_obs)

demand = vec(remove_missing(xf_demand["Sheet1"][sprintf1("B2:O%d", num_t+1)])); # [z+t]
lambda_obs = vec(remove_missing(xf_da_prices["Sheet1"][sprintf1("B2:O%d", num_t+1)])); # [z+t]
np_obs = vec(remove_missing(xf_netpos["Sheet1"][sprintf1("B2:O%d", num_t+1)])); # [z+t]

demand = convert(Vector{Float64}, demand)
lambda_obs = convert(Vector{Float64}, lambda_obs)
np_obs = convert(Vector{Float64}, np_obs)

ptdf_z = vec(remove_missing(xf_ptdf["Sheet1"][sprintf1("D2:Q%d", num_t*num_j+1)]))
ram = vec(remove_missing(xf_ptdf["Sheet1"][sprintf1("C2:C%d", num_t*num_j+1)]))

ptdf_z = convert(Vector{Float64}, ptdf_z)
ram = convert(Vector{Float64}, ram)

g_max = vec(xf_capacities["Sheet1"]["B2:O11"]); # [z+tech]
g_max_t = zeros(num_z*num_tech*num_t); # [z+tech+t]
for z in 1:num_z
    for t in 1:num_t
        for tech in 1:num_tech
            g_max_t[num_t*num_tech*(z-1)+num_t*(tech-1)+t] = g_max[num_tech*(z-1)+tech] 
        end
    end
end

ren_gen = vec(remove_missing(xf_ren_gen["Sheet1"][sprintf1("B2:O%d", num_t+1)])); # [z+t]
ren_gen = convert(Vector{Float64}, ren_gen);