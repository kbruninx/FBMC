using NPZ
using SparseArrays
using Plots
using DataFrames, XLSX
using Dates
using JuMP, HiGHS, Ipopt, Alpine, Juniper
using JLD

zones =  ["ALBE", "ALDE", "AT", "BE", "CZ", "DE_LU", "FR", "HR", "HU", "NL", "PL", "RO", "SI", "SK"]
countries =  ["AT", "BE", "CZ", "DE", "FR", "HR", "HU", "NL", "PL", "RO", "SI", "SK"]

df_ptdf = DataFrame(XLSX.readtable("./flow_based_domain/ptdf_z_obs.xlsx", "Sheet1"))
df_ptdf.DateTime = DateTime.(df_ptdf.DateTime)
df_ptdf = df_ptdf[df_ptdf.DateTime .<= DateTime(2023, 3, 30, 23), :]

df_timestamps = DataFrame(XLSX.readtable("./flow_based_domain/timestamps.xlsx", "Sheet1"))
df_timestamps.DateTime = DateTime.(df_timestamps.DateTime)

df_contingencies = DataFrame(XLSX.readtable("./flow_based_domain/contingencies.xlsx", "Sheet1"))
df_cnecs = DataFrame(XLSX.readtable("./flow_based_domain/ptdf_z_obs_means.xlsx", "Sheet1"))
df_line_edge_map = DataFrame(XLSX.readtable("./flow_based_domain/line_edge_map.xlsx", "Sheet1"))
df_grid = DataFrame(XLSX.readtable("./flow_based_domain/grid.xlsx", "Sheet1"))

H_mat = npzread("./flow_based_domain/H.npy")
L_mat = npzread("./flow_based_domain/L.npy")

GSK1_P = npzread("./flow_based_domain/gsk1_matrix_p.npy")
GSK2_P = npzread("./flow_based_domain/gsk2_matrix_p.npy")
GSK3_P = npzread("./flow_based_domain/gsk3_matrix_p.npy")

#start_date = DateTime(2022, 9, 11)
start_date = DateTime(2023, 2, 1)
end_date = DateTime(2023, 3, 1)
start_t = findfirst(==(start_date), df_timestamps.DateTime)
start_t_bc = start_t - Dates.value(convert(Dates.Hour, DateTime(2022, 9, 11) - DateTime(2022, 9, 1)))
hour_count = Dates.value(convert(Dates.Hour, end_date - start_date))

GSK1_P = GSK1_P[start_t:(start_t+hour_count-1), :, :]
GSK2_P = GSK2_P[start_t:(start_t+hour_count-1), :, :]
GSK3_P = GSK3_P[start_t:(start_t+hour_count-1), :, :]

GSK4_P = load("./flow_based_domain/gsk4_matrix_p.jld")["data"]
GSK5_P = load("./flow_based_domain/gsk5_matrix_p.jld")["data"]

GSK4_P = GSK4_P[start_t_bc:(start_t_bc+hour_count-1), :, :]
GSK5_P = GSK5_P[start_t_bc:(start_t_bc+hour_count-1), :, :]

O = npzread("./flow_based_domain/omega_matrix.npy")

PTDF_N_C = Array{Matrix{Float64}}(undef, size(df_contingencies)[1])
for i = 1:size(df_contingencies)[1]
    PTDF_N_C[i] = npzread("./flow_based_domain/ptdf_n/ptdf_n_$(i-1).npy")
end

num_l = size(H_mat)[1]
num_gsk_strategy = 5