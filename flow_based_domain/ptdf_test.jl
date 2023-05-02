using NPZ
using SparseArrays
using Plots
using DataFrames, XLSX
using Dates

zones =  ["ALBE", "ALDE", "AT", "BE", "CZ", "DE_LU", "FR", "HR", "HU", "NL", "PL", "RO", "SI", "SK"]

df_ptdf = DataFrame(XLSX.readtable("./flow_based_domain/ptdf_z_obs.xlsx", "Sheet1"))
df_ptdf.DateTime = DateTime.(df_ptdf.DateTime)

GSK1_P = npzread("./flow_based_domain/gsk1_matrix_p.npy")
GSK2_P = npzread("./flow_based_domain/gsk2_matrix_p.npy")
GSK3_P = npzread("./flow_based_domain/gsk3_matrix_p.npy")
GSK5_P = npzread("./flow_based_domain/gsk5_matrix_p.npy")
O = npzread("./flow_based_domain/omega_matrix.npy")

M = zeros(size(O)[2], size(O)[1])
X = zeros(size(O)[3], size(GSK1_P)[2])
X[1, :] .= 1
#X[1, :] .= 0.6
#X[2, :] .= 0.2
#X[3, :] .= 0.2
#X[:, :] .= 0.2

for p in 1:size(GSK1_P)[2]
    M[:, p] = O[p, :, :] * X[:, p]
end

"""
PTDF1_Z_T = Array{Matrix{Float64}}(undef, size(GSK1_P)[1])
for t = 1:size(GSK1_P)[1]
    PTDF1_Z_T[t] = PTDF_N * M * GSK1_P[t, :, :]
    println(t)
end

PTDF2_Z_T = Array{Matrix{Float64}}(undef, size(GSK2_P)[1])
for t = 1:size(GSK2_P)[1]
    PTDF2_Z_T[t] = PTDF_N * M * GSK2_P[t, :, :]
    println(t)
end

PTDF3_Z_T = Array{Matrix{Float64}}(undef, size(GSK3_P)[1])
for t = 1:size(GSK3_P)[1]
    PTDF3_Z_T[t] = PTDF_N * M * GSK3_P[t, :, :]
    println(t)
end
"""

"""
function plot_ptdf(line_id, zone, centring_value)
    ptdf1_line_t = Array{Float64}(undef, size(GSK1_P)[1])
    for t = 1:size(GSK1_P)[1]
        ptdf1_line_t[t] = PTDF1_Z_T[t][line_id, zone] + centring_value
    end

    ptdf2_line_t = Array{Float64}(undef, size(GSK2_P)[1])
    for t = 1:size(GSK2_P)[1]
        ptdf2_line_t[t] = PTDF2_Z_T[t][line_id, zone] + centring_value
    end

    ptdf3_line_t = Array{Float64}(undef, size(GSK3_P)[1])
    for t = 1:size(GSK3_P)[1]
        ptdf3_line_t[t] = PTDF3_Z_T[t][line_id, zone] + centring_value
    end
    plot!(twiny(), [ptdf1_line_t, ptdf2_line_t, ptdf3_line_t])
end

plot(df_ptdf[df_ptdf.line_id .== 815, :DateTime], df_ptdf[df_ptdf.line_id .== 815, :HU])
plot_ptdf(800, 9, -0.026198)
"""

function generate_PTDF(line, edge_number, contingency, zone, centring_value)
    ptdf_z_obs_x = df_ptdf[(df_ptdf.line_id .== line) .& (df_ptdf.contingency .== contingency), :DateTime]
    ptdf_z_obs_y = df_ptdf[(df_ptdf.line_id .== line) .& (df_ptdf.contingency .== contingency), zone]
    plot(ptdf_z_obs_x, ptdf_z_obs_y)

    PTDF_N = npzread("./flow_based_domain/ptdf_n/ptdf_n_$contingency.npy")

    PTDF1_Z_T = Array{Matrix{Float64}}(undef, size(GSK1_P)[1])
    for t = 1:size(GSK1_P)[1]
        PTDF1_Z_T[t] = PTDF_N * M * GSK1_P[t, :, :]
        println(t)
    end

    PTDF2_Z_T = Array{Matrix{Float64}}(undef, size(GSK2_P)[1])
    for t = 1:size(GSK2_P)[1]
        PTDF2_Z_T[t] = PTDF_N * M * GSK2_P[t, :, :]
        println(t)
    end

    PTDF3_Z_T = Array{Matrix{Float64}}(undef, size(GSK3_P)[1])
    for t = 1:size(GSK3_P)[1]
        PTDF3_Z_T[t] = PTDF_N * M * GSK3_P[t, :, :]
        println(t)
    end

    PTDF5_Z_T = Array{Matrix{Float64}}(undef, size(GSK5_P)[1])
    for t = 1:size(GSK5_P)[1]
        PTDF5_Z_T[t] = PTDF_N * M * GSK5_P[t, :, :]
        println(t)
    end

    zone_i = findall(zones.==zone)[1]
    ptdf1_line_t = Array{Float64}(undef, size(GSK1_P)[1])
    for t = 1:size(GSK1_P)[1]
        ptdf1_line_t[t] = PTDF1_Z_T[t][edge_number + 1, zone_i]# + centring_value
    end

    ptdf2_line_t = Array{Float64}(undef, size(GSK2_P)[1])
    for t = 1:size(GSK2_P)[1]
        ptdf2_line_t[t] = PTDF2_Z_T[t][edge_number + 1, zone_i]# + centring_value
    end

    ptdf3_line_t = Array{Float64}(undef, size(GSK3_P)[1])
    for t = 1:size(GSK3_P)[1]
        ptdf3_line_t[t] = PTDF3_Z_T[t][edge_number + 1, zone_i]# + centring_value
    end

    ptdf5_line_t = Array{Float64}(undef, size(GSK5_P)[1])
    for t = 1:size(GSK5_P)[1]
        ptdf5_line_t[t] = PTDF5_Z_T[t][edge_number + 1, zone_i]# + centring_value
    end

    plot!(twiny(), [ptdf1_line_t, ptdf2_line_t, ptdf3_line_t, ptdf5_line_t])
    #plot([ptdf1_line_t, ptdf2_line_t, ptdf3_line_t, ptdf5_line_t])
end