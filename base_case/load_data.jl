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