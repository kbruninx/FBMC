using Dates
using JuMP, HiGHS
using DataFrames, XLSX
using LinearAlgebra
using Statistics
using SparseArrays
using Formatting
using Polyhedra, CDDLib
using JLD
using LaTeXStrings
using Plots

zones =  ["AT", "BE", "CZ", "DE_LU", "FR", "HR", "HU", "NL", "PL", "RO", "SI", "SK"]

df_ptdf_obs = DataFrame(XLSX.readtable("./flow_based_domain/ptdf_z_obs.xlsx", "Sheet1"))
df_ptdf_obs.DateTime = DateTime.(df_ptdf_obs.DateTime)

df_ptdf_10 = DataFrame(XLSX.readtable("./flow_based_domain/ptdf_z_calc_nuts.xlsx", "Sheet1"))
df_ptdf_10.DateTime = DateTime.(df_ptdf_10.DateTime)

df_ptdf_5 = DataFrame(XLSX.readtable("./flow_based_domain/ptdf_z_calc_nuts_feb_5.xlsx", "Sheet1"))
df_ptdf_5.DateTime = DateTime.(df_ptdf_5.DateTime)

df_ptdf_n = DataFrame(XLSX.readtable("./flow_based_domain/ptdf_z_naive.xlsx", "Sheet1"))
df_ptdf_n.DateTime = DateTime.(df_ptdf_n.DateTime)

df_timestamps = DataFrame(XLSX.readtable("./flow_based_domain/timestamps.xlsx", "Sheet1"))
df_timestamps.DateTime = DateTime.(df_timestamps.DateTime)

df_np_obs = DataFrame(XLSX.readtable("./cost_curves/data/net_positions.xlsx", "Sheet1"))
df_np_obs.DateTime = DateTime.(df_np_obs.DateTime)

start_date = DateTime(2023, 2, 1)
end_date = DateTime(2023, 3, 1)
df_timestamps = df_timestamps[(df_timestamps.DateTime .>= start_date) .& (df_timestamps.DateTime .<= end_date), :DateTime]
