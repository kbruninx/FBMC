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

start_date = DateTime(2022, 9, 11)
end_date = DateTime(2023, 3, 1)
start_t = findfirst(==(start_date), df_timestamps.DateTime)
hour_count = Dates.value(convert(Dates.Hour, end_date - start_date))

GSK1_P = GSK1_P[start_t:(start_t+hour_count-1), :, :]
GSK2_P = GSK2_P[start_t:(start_t+hour_count-1), :, :]
GSK3_P = GSK3_P[start_t:(start_t+hour_count-1), :, :]

GSK4_P = load("./flow_based_domain/gsk4_matrix_p.jld")["data"]
GSK5_P = load("./flow_based_domain/gsk5_matrix_p.jld")["data"]

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

#function generate_PTDF(line, edge_number, contingency, zone, centring_value)
    line = 2813
    edge_number = 792
    contingency = 16
    zone = "SK"
    centring_value = 0
    ptdf_z_obs_x = df_ptdf[(df_ptdf.line_id .== line) .& (df_ptdf.contingency .== contingency) .& (df_ptdf.DateTime .>= start_date) .& (df_ptdf.DateTime .<= end_date), :DateTime]
    ptdf_z_obs_y = df_ptdf[(df_ptdf.line_id .== line) .& (df_ptdf.contingency .== contingency) .& (df_ptdf.DateTime .>= start_date) .& (df_ptdf.DateTime .<= end_date), zone]
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

    PTDF4_Z_T = Array{Matrix{Float64}}(undef, size(GSK4_P)[1])
    for t = 1:size(GSK5_P)[1]
        PTDF4_Z_T[t] = PTDF_N * M * GSK4_P[t, :, :]
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

    ptdf4_line_t = Array{Float64}(undef, size(GSK4_P)[1])
    for t = 1:size(GSK4_P)[1]
        ptdf4_line_t[t] = PTDF4_Z_T[t][edge_number + 1, zone_i]# + centring_value
    end

    ptdf5_line_t = Array{Float64}(undef, size(GSK5_P)[1])
    for t = 1:size(GSK5_P)[1]
        ptdf5_line_t[t] = PTDF5_Z_T[t][edge_number + 1, zone_i]# + centring_value
    end

    plot!(twiny(), [ptdf1_line_t, ptdf2_line_t, ptdf3_line_t, ptdf4_line_t, ptdf5_line_t])
    #plot([ptdf1_line_t, ptdf2_line_t, ptdf3_line_t, ptdf5_line_t])
#end