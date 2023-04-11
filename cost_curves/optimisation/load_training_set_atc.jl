using JuMP, HiGHS
#using Gurobi
using DataFrames, XLSX
using LinearAlgebra
using Statistics
using SparseArrays
using Formatting

function remove_missing(mx)
    mx[ismissing.(mx)] .= 0.0
    return mx
end

fbmc_zones =  ["ALBE", "ALDE", "AT", "BE", "CZ", "DE_LU", "FR", "HR", "HU", "NL", "PL", "RO", "SI", "SK"]
non_fbmc_zones =  ["CH", "GB", "ES", "IT_NORD"] # ["CH", "GB", "ES", "IT_NORD", "DK_1", "NO_2"]
all_zones = vcat(fbmc_zones, non_fbmc_zones)

borders = [
    ["FR", "CH"],
    ["FR", "ES"],
    ["FR", "GB"],
    ["FR", "IT_NORD"],
    ["AT", "IT_NORD"],
    ["AT", "CH"],
    ["BE", "GB"],
    ["NL", "GB"],
    ["SI", "IT_NORD"],
    #["NL", "DK_1"],
    #["NL", "NO_2"],
]

plant_types = [
    "biomass", "brown_coal", "coal_gas", "natural_gas", "hard_coal", "oil", "hydro", 
    "nuclear", "waste", "other"
] 

df_atc = DataFrame(XLSX.readtable("./data/training/atc.xlsx", "Sheet1"))

xf_da_prices = XLSX.readxlsx("./data/training/day_ahead_prices.xlsx")
xf_demand = XLSX.readxlsx("./data/training/demand.xlsx")
xf_fuel_prices = XLSX.readxlsx("./data/training/fuel_prices.xlsx")
xf_generation = XLSX.readxlsx("./data/training/generation.xlsx")
xf_ren_gen = XLSX.readxlsx("./data/training/renewable_generation.xlsx")
xf_netpos = XLSX.readxlsx("./data/training/net_positions.xlsx")
xf_ptdf = XLSX.readxlsx("./data/training/ptdfs.xlsx")
xf_capacities = XLSX.readxlsx("./data/installed_capacities_corrected.xlsx")

xf_da_prices_non_fbmc = XLSX.readxlsx("./data/training/day_ahead_prices_non_fbmc.xlsx")
xf_demand_non_fbmc = XLSX.readxlsx("./data/training/demand_non_fbmc.xlsx")
xf_generation_non_fbmc = XLSX.readxlsx("./data/training/generation_non_fbmc.xlsx")
xf_ren_gen_non_fbmc = XLSX.readxlsx("./data/training/renewable_generation_non_fbmc.xlsx")
xf_capacities_non_fbmc = XLSX.readxlsx("./data/installed_capacities_non_fbmc.xlsx")


num_z = 14 # including ALBE and ALDE
num_z_non_fbmc = 4 #6
num_atc_border = 9 #11

num_train_t = 2906 # size of training set
num_tech = 10
num_j = 133 # maximum amont on CNEs at a given time in the training dataset

coal_prices_g = vec(remove_missing(xf_fuel_prices["Sheet1"][sprintf1("B2:B%d", num_train_t+1)]))
oil_prices_g = vec(remove_missing(xf_fuel_prices["Sheet1"][sprintf1("C2:C%d", num_train_t+1)]))
gas_prices_g = vec(remove_missing(xf_fuel_prices["Sheet1"][sprintf1("D2:D%d", num_train_t+1)]))
eua_prices_g = vec(remove_missing(xf_fuel_prices["Sheet1"][sprintf1("E2:E%d", num_train_t+1)]))

coal_prices_g = convert(Vector{Float64}, coal_prices_g)
oil_prices_g = convert(Vector{Float64}, oil_prices_g)
gas_prices_g = convert(Vector{Float64}, gas_prices_g)
eua_prices_g = convert(Vector{Float64}, eua_prices_g)

at_obs_g = remove_missing(xf_generation["AT"][sprintf1("B2:K%d", num_train_t+1)])
be_obs_g = remove_missing(xf_generation["BE"][sprintf1("B2:K%d", num_train_t+1)])
cz_obs_g = remove_missing(xf_generation["CZ"][sprintf1("B2:K%d", num_train_t+1)])
de_obs_g = remove_missing(xf_generation["DE_LU"][sprintf1("B2:K%d", num_train_t+1)])
fr_obs_g = remove_missing(xf_generation["FR"][sprintf1("B2:K%d", num_train_t+1)])
hr_obs_g = remove_missing(xf_generation["HR"][sprintf1("B2:K%d", num_train_t+1)])
hu_obs_g = remove_missing(xf_generation["HU"][sprintf1("B2:K%d", num_train_t+1)])
nl_obs_g = remove_missing(xf_generation["NL"][sprintf1("B2:K%d", num_train_t+1)])
pl_obs_g = remove_missing(xf_generation["PL"][sprintf1("B2:K%d", num_train_t+1)])
ro_obs_g = remove_missing(xf_generation["RO"][sprintf1("B2:K%d", num_train_t+1)])
si_obs_g = remove_missing(xf_generation["SI"][sprintf1("B2:K%d", num_train_t+1)])
sk_obs_g = remove_missing(xf_generation["SK"][sprintf1("B2:K%d", num_train_t+1)])

ch_obs_g = remove_missing(xf_generation_non_fbmc["CH"][sprintf1("B2:K%d", num_train_t+1)])
gb_obs_g = remove_missing(xf_generation_non_fbmc["GB"][sprintf1("B2:K%d", num_train_t+1)])
es_obs_g = remove_missing(xf_generation_non_fbmc["ES"][sprintf1("B2:K%d", num_train_t+1)])
it_obs_g = remove_missing(xf_generation_non_fbmc["IT_NORD"][sprintf1("B2:K%d", num_train_t+1)])
dk_obs_g = remove_missing(xf_generation_non_fbmc["DK_1"][sprintf1("B2:K%d", num_train_t+1)])
no_obs_g = remove_missing(xf_generation_non_fbmc["NO_2"][sprintf1("B2:K%d", num_train_t+1)])

demand_g = remove_missing(xf_demand["Sheet1"][sprintf1("B2:O%d", num_train_t+1)]) # [z+t]
#demand_g_non_fbmc = remove_missing(xf_demand_non_fbmc["Sheet1"][sprintf1("B2:G%d", num_train_t+1)]) # [z+t]
demand_g_non_fbmc = remove_missing(xf_demand_non_fbmc["Sheet1"][sprintf1("B2:E%d", num_train_t+1)]) 

lambda_obs_g = remove_missing(xf_da_prices["Sheet1"][sprintf1("B2:O%d", num_train_t+1)]) # [z+t]
#lambda_obs_g_non_fbmc = remove_missing(xf_da_prices_non_fbmc["Sheet1"][sprintf1("B2:G%d", num_train_t+1)]) # [z+t]
lambda_obs_g_non_fbmc = remove_missing(xf_da_prices_non_fbmc["Sheet1"][sprintf1("B2:E%d", num_train_t+1)]) # [z+t]

np_obs_g = remove_missing(xf_netpos["Sheet1"][sprintf1("B2:O%d", num_train_t+1)]) # [z+t]

ptdf_z_g = remove_missing(xf_ptdf["Sheet1"][sprintf1("D2:Q%d", num_train_t*num_j+1)])
ram_g = remove_missing(xf_ptdf["Sheet1"][sprintf1("C2:C%d", num_train_t*num_j+1)])

g_max_g = vec(xf_capacities["Sheet1"]["B2:O11"]) # [z+tech]
#g_max_g_non_fbmc = vec(xf_capacities_non_fbmc["Sheet1"]["B2:G11"]) # [z+tech]
g_max_g_non_fbmc = vec(xf_capacities_non_fbmc["Sheet1"]["B2:E11"]) # [z+tech]

ren_gen_g = remove_missing(xf_ren_gen["Sheet1"][sprintf1("B2:O%d", num_train_t+1)]) # [z+t]
#ren_gen_g_non_fbmc = remove_missing(xf_ren_gen_non_fbmc["Sheet1"][sprintf1("B2:G%d", num_train_t+1)]) # [z+t]
ren_gen_g_non_fbmc = remove_missing(xf_ren_gen_non_fbmc["Sheet1"][sprintf1("B2:E%d", num_train_t+1)]) # [z+t]