alpha = 0.04
sigma = load("./flow_based_domain/sigma_nuts.jld")["data"]

M = zeros(size(O)[2], size(O)[1])
X = zeros(size(O)[3], size(GSK1_P)[2])
X[1, :] .= 1

for p in 1:size(GSK1_P)[2]
    M[:, p] = O[p, :, :] * X[:, p]
end

df_ptdf_calc = DataFrame(vcat([[],[],[],[],[],[]], [[] for z in zones]), vcat(["DateTime", "line_id", "contingency", "fref", "frm", "ram"], zones))

for t in 1:hour_count
#for t in 1:1
    println(t)

    PTDF_Z_s = [
        PTDF_N_C[1] * M * GSK1_P[t, :, :];;;
        PTDF_N_C[1] * M * GSK2_P[t, :, :];;;
        PTDF_N_C[1] * M * GSK3_P[t, :, :];;;
        PTDF_N_C[1] * M * GSK4_P[t, :, :];;;
        PTDF_N_C[1] * M * GSK5_P[t, :, :];;;
    ]

    PTDF_Z_z = Array{Matrix, 1}(undef, size(countries)[1])
    for z in 1:size(countries)[1]
        PTDF_Z_z[z] = sigma[1, z] * PTDF_Z_s[:, :, 1]

        for s in 2:size(sigma[:, z])[1]
            PTDF_Z_z[z] += sigma[s, z] * PTDF_Z_s[:, :, s]
        end
    end

    for l = 1:num_l
        line_id = df_line_edge_map[df_line_edge_map.edge .== l - 1, :line_id][1]
        zone = df_grid[df_grid.line_id .== line_id, :zone][1]
        #tieline = df_grid[df_grid.line_id .== line_id, :tieline][1] && df_grid[df_grid.line_id .== line_id, :eic][1][1:4] == "10T-"

        z = findfirst(==(zone), countries)
        PTDF_Z = PTDF_Z_z[z]

        ptdf_l_z = PTDF_Z[l, z]
        z_others = 1:size(countries)[1]
        z_others = z_others[z_others .!= z]
        ptdf_l_z_others = [PTDF_Z[l, z_i] for z_i in z_others]

        is_cne = false
        #if tieline
        #    is_cne = true
        #else
        for ptdf_l_z_i in ptdf_l_z_others
            if abs(ptdf_l_z - ptdf_l_z_i) >= alpha
                is_cne = true
            end
        end
        #end

        if is_cne
            timestamp = df_timestamps[start_t + t, :DateTime]
            push!(df_ptdf_calc, vcat([timestamp, line_id, 0, 0, 0, 0], PTDF_Z[l, :]))
        end
    end
end

XLSX.writetable("./flow_based_domain/ptdf_z_calc_nuts_feb_4.xlsx", "Sheet1" => df_ptdf_calc)