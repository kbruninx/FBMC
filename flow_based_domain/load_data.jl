using NPZ
using SparseArrays
using Plots
using DataFrames, XLSX
using Dates
using JuMP, HiGHS, Ipopt, Alpine

zones =  ["ALBE", "ALDE", "AT", "BE", "CZ", "DE_LU", "FR", "HR", "HU", "NL", "PL", "RO", "SI", "SK"]

df_ptdf = DataFrame(XLSX.readtable("./flow_based_domain/ptdf_z_obs.xlsx", "Sheet1"))
df_ptdf.DateTime = DateTime.(df_ptdf.DateTime)

df_timestamps = DataFrame(XLSX.readtable("./flow_based_domain/timestamps.xlsx", "Sheet1"))
df_timestamps.DateTime = DateTime.(df_timestamps.DateTime)

df_contingencies = DataFrame(XLSX.readtable("./flow_based_domain/contingencies.xlsx", "Sheet1"))
df_cnecs = DataFrame(XLSX.readtable("./flow_based_domain/ptdf_z_obs_means.xlsx", "Sheet1"))
df_line_edge_map = DataFrame(XLSX.readtable("./flow_based_domain/line_edge_map.xlsx", "Sheet1"))

GSK1_P = npzread("./flow_based_domain/gsk1_matrix_p.npy")
GSK2_P = npzread("./flow_based_domain/gsk2_matrix_p.npy")
GSK3_P = npzread("./flow_based_domain/gsk3_matrix_p.npy")
O = npzread("./flow_based_domain/omega_matrix.npy")

PTDF_N_C = Array{Matrix{Float64}}(undef, size(df_contingencies)[1])
for i = 1:size(df_contingencies)[1]
    PTDF_N_C[i] = npzread("./flow_based_domain/ptdf_n/ptdf_n_$(i-1).npy")
end

num_gsk_strategy = 3