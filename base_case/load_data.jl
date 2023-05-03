using NPZ
using Dates
using JuMP, HiGHS
using DataFrames, XLSX
using LinearAlgebra
using Statistics
using SparseArrays
using Formatting

function remove_missing(mx)
    mx[ismissing.(mx)] .= 0.0
    return mx
end


non_fbmc_zones =  ["CH", "GB", "ES", "IT_NORD"]
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
    ["DE_LU", "CH"],
    ["IT_NORD", "CH"],
]

plant_types = [
    "biomass", "brown_coal", "coal_gas", "natural_gas", "hard_coal", "oil", "hydro", 
    "nuclear", "waste", "other"
] 

df_atc = DataFrame(XLSX.readtable("./data/atc.xlsx", "Sheet1"))

xf_fuel_prices = XLSX.readxlsx("./data/fuel_prices.xlsx")

xf_demand = XLSX.readxlsx("./cost_curves/data/demand.xlsx")
xf_ren_gen = XLSX.readxlsx("./data/renewable_generation.xlsx")
xf_capacities = XLSX.readxlsx("./data/installed_capacities.xlsx")
xf_generation_outages = XLSX.readxlsx("./data/generation_outages.xlsx")

xf_demand_non_fbmc = XLSX.readxlsx("./data/demand_non_fbmc.xlsx")
xf_ren_gen_non_fbmc = XLSX.readxlsx("./data/renewable_generation_non_fbmc.xlsx")
xf_capacities_non_fbmc = XLSX.readxlsx("./data/installed_capacities_non_fbmc_corrected.xlsx")
xf_generation_outages_non_fbmc = XLSX.readxlsx("./data/generation_outages_non_fbmc.xlsx")


num_z = 14 # including ALBE and ALDE
num_z_non_fbmc = 4
num_atc_border = 11

num_train_t = 5064 # size of training set
num_tech = 10
num_j = 133 # maximum amount of CNEs at a given time in the training dataset

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

at_gen_out_g = remove_missing(xf_generation_outages["AT"][sprintf1("B2:K%d", num_train_t+1)])
be_gen_out_g = remove_missing(xf_generation_outages["BE"][sprintf1("B2:K%d", num_train_t+1)])
cz_gen_out_g = remove_missing(xf_generation_outages["CZ"][sprintf1("B2:K%d", num_train_t+1)])
de_gen_out_g = remove_missing(xf_generation_outages["DE_LU"][sprintf1("B2:K%d", num_train_t+1)])
fr_gen_out_g = remove_missing(xf_generation_outages["FR"][sprintf1("B2:K%d", num_train_t+1)])
hr_gen_out_g = remove_missing(xf_generation_outages["HR"][sprintf1("B2:K%d", num_train_t+1)])
hu_gen_out_g = remove_missing(xf_generation_outages["HU"][sprintf1("B2:K%d", num_train_t+1)])
nl_gen_out_g = remove_missing(xf_generation_outages["NL"][sprintf1("B2:K%d", num_train_t+1)])
pl_gen_out_g = remove_missing(xf_generation_outages["PL"][sprintf1("B2:K%d", num_train_t+1)])
ro_gen_out_g = remove_missing(xf_generation_outages["RO"][sprintf1("B2:K%d", num_train_t+1)])
si_gen_out_g = remove_missing(xf_generation_outages["SI"][sprintf1("B2:K%d", num_train_t+1)])
sk_gen_out_g = remove_missing(xf_generation_outages["SK"][sprintf1("B2:K%d", num_train_t+1)])
ch_gen_out_g = remove_missing(xf_generation_outages_non_fbmc["CH"][sprintf1("B2:K%d", num_train_t+1)])
gb_gen_out_g = remove_missing(xf_generation_outages_non_fbmc["GB"][sprintf1("B2:K%d", num_train_t+1)])
es_gen_out_g = remove_missing(xf_generation_outages_non_fbmc["ES"][sprintf1("B2:K%d", num_train_t+1)])
it_gen_out_g = remove_missing(xf_generation_outages_non_fbmc["IT_NORD"][sprintf1("B2:K%d", num_train_t+1)])

demand_g = remove_missing(xf_demand["Sheet1"][sprintf1("B2:O%d", num_train_t+1)]) # [z+t]
demand_g_non_fbmc = remove_missing(xf_demand_non_fbmc["Sheet1"][sprintf1("B2:E%d", num_train_t+1)]) 

lambda_obs_g = remove_missing(xf_da_prices["Sheet1"][sprintf1("B2:O%d", num_train_t+1)]) # [z+t]
lambda_obs_g_non_fbmc = remove_missing(xf_da_prices_non_fbmc["Sheet1"][sprintf1("B2:E%d", num_train_t+1)]) # [z+t]

np_obs_g = remove_missing(xf_netpos["Sheet1"][sprintf1("B2:O%d", num_train_t+1)]) # [z+t]

ptdf_z_g = remove_missing(xf_ptdf["Sheet1"][sprintf1("D2:Q%d", num_train_t*num_j+1)])
ram_g = remove_missing(xf_ptdf["Sheet1"][sprintf1("C2:C%d", num_train_t*num_j+1)])

g_max_g = vec(xf_capacities["Sheet1"]["B2:O11"]) # [z+tech]
g_max_g_non_fbmc = vec(xf_capacities_non_fbmc["Sheet1"]["B2:E11"]) # [z+tech]

ren_gen_g = remove_missing(xf_ren_gen["Sheet1"][sprintf1("B2:O%d", num_train_t+1)]) # [z+t]
ren_gen_g_non_fbmc = remove_missing(xf_ren_gen_non_fbmc["Sheet1"][sprintf1("B2:E%d", num_train_t+1)]) # [z+t]

timestamps = DateTime.(vec(remove_missing(xf_demand["Sheet1"][sprintf1("A2:A%d", num_train_t+1)]))) 


# USEFUL PART

zones =  ["AT", "BE", "CZ", "DE", "FR", "HR", "HU", "NL", "PL", "RO", "SI", "SK"]

df_timestamps = DataFrame(XLSX.readtable("./flow_based_domain/timestamps.xlsx", "Sheet1"))
df_timestamps.DateTime = DateTime.(df_timestamps.DateTime)

df_plants = DataFrame(XLSX.readtable("./flow_based_domain/plants.xlsx", "Sheet1"))
df_substations = DataFrame(XLSX.readtable("./flow_based_domain/substations.xlsx", "Sheet1"))
df_grid = DataFrame(XLSX.readtable("./flow_based_domain/grid.xlsx", "Sheet1"))

H_mat = npzread("./flow_based_domain/H.npy")
L_mat = npzread("./flow_based_domain/L.npy")

df_demand = DataFrame(XLSX.readtable("./cost_curves/data/demand.xlsx", "Sheet1"))
demand_t_z = remove_missing(Matrix(df_demand[:, 3:14]))

df_ren_gen = DataFrame(XLSX.readtable("./cost_curves/data/renewable_generation.xlsx", "Sheet1"))
ren_gen_t_z = remove_missing(Matrix(df_ren_gen[:, 3:14]))

df_line_edge_map = DataFrame(XLSX.readtable("./flow_based_domain/line_edge_map.xlsx", "Sheet1"))
df_substation_node_map = DataFrame(XLSX.readtable("./flow_based_domain/substation_node_map.xlsx", "Sheet1"))