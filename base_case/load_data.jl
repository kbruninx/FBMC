using NPZ
using Dates
using JuMP, HiGHS, Juniper, Ipopt
using DataFrames, XLSX
using LinearAlgebra
using Statistics
using SparseArrays
using Formatting
using JLD

num_train_t = 7948 # size of training set

function remove_missing(mx)
    mx[ismissing.(mx)] .= 0.0
    return mx
end

zones =  ["AT", "BE", "CZ", "DE", "FR", "HR", "HU", "NL", "PL", "RO", "SI", "SK"]

plant_types = [
    "biomass", "brown_coal", "coal_gas", "natural_gas", "hard_coal", "oil", "hydro", 
    "nuclear", "waste", "other"
] 

df_timestamps = DataFrame(XLSX.readtable("./flow_based_domain/timestamps_july.xlsx", "Sheet1"))
df_timestamps.DateTime = DateTime.(df_timestamps.DateTime)

df_plants = DataFrame(XLSX.readtable("./flow_based_domain/plants.xlsx", "Sheet1"))
df_substations = DataFrame(XLSX.readtable("./flow_based_domain/substations.xlsx", "Sheet1"))
df_grid = DataFrame(XLSX.readtable("./flow_based_domain/grid.xlsx", "Sheet1"))

plant_weights = 1 .- normalize(df_plants.construction_year)

H_mat = npzread("./flow_based_domain/H.npy")
L_mat = npzread("./flow_based_domain/L.npy")
availability_matrix = npzread("./flow_based_domain/availability_matrix_july.npy")

df_demand = DataFrame(XLSX.readtable("./cost_curves/data-july/demand.xlsx", "Sheet1"))
demand_t_z = remove_missing(Matrix(df_demand[:, 3:14]))

df_ren_gen = DataFrame(XLSX.readtable("./cost_curves/data-july/renewable_generation.xlsx", "Sheet1"))
ren_gen_t_z = remove_missing(Matrix(df_ren_gen[:, 3:14]))

df_line_edge_map = DataFrame(XLSX.readtable("./flow_based_domain/line_edge_map.xlsx", "Sheet1"))
df_substation_node_map = DataFrame(XLSX.readtable("./flow_based_domain/substation_node_map.xlsx", "Sheet1"))

xf_fuel_prices = XLSX.readxlsx("./cost_curves/data-july/fuel_prices.xlsx")

coal_prices_g = vec(remove_missing(xf_fuel_prices["Sheet1"][sprintf1("B2:B%d", num_train_t+1)]))
oil_prices_g = vec(remove_missing(xf_fuel_prices["Sheet1"][sprintf1("C2:C%d", num_train_t+1)]))
gas_prices_g = vec(remove_missing(xf_fuel_prices["Sheet1"][sprintf1("D2:D%d", num_train_t+1)]))
eua_prices_g = vec(remove_missing(xf_fuel_prices["Sheet1"][sprintf1("E2:E%d", num_train_t+1)]))

coal_prices_g = convert(Vector{Float64}, coal_prices_g)
oil_prices_g = convert(Vector{Float64}, oil_prices_g)
gas_prices_g = convert(Vector{Float64}, gas_prices_g)
eua_prices_g = convert(Vector{Float64}, eua_prices_g)

#ref_flow_obs = npzread("./base_case/ref_flow_obs.npy")

num_p = size(df_plants)[1]
num_n = size(H_mat)[2]
num_l = size(H_mat)[1]
num_z = 12
num_t = 24 # running base case on a daily basis
num_tech = 10