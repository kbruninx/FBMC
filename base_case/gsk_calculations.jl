fbmc_zones =  ["ALBE", "ALDE", "AT", "BE", "CZ", "DE_LU", "FR", "HR", "HU", "NL", "PL", "RO", "SI", "SK"]

ref_g = load(string("./base_case/ref_g.jld"))["data"]
ref_g_nuts = load(string("./base_case/ref_g_nuts.jld"))["data"]

start_date = DateTime(2022, 9, 11)
end_date = DateTime(2023, 3, 1)
hour_count = Dates.value(convert(Dates.Hour, end_date - start_date))

# GSK4: generation level

# denominator
gsk4_p_plants_z = Dict()
for zone in fbmc_zones
    plants_z = findall(==(zone), df_plants.zone)
    gsk4_p_plants_z[zone] = zeros(hour_count)
    for t in 1:hour_count
        sum_plants_z = 0
        for p in plants_z
            sum_plants_z += ref_g[t, p]
        end
        gsk4_p_plants_z[zone][t] = sum_plants_z
    end
end

# nominator
gsk4_matrix = zeros(hour_count, num_p, size(fbmc_zones)[1])
for p in 1:num_p
    zone = df_plants[p, "zone"]
    for t in 1:hour_count
        if gsk4_p_plants_z[zone][t] > 0
            gsk4_matrix[t, p, findfirst(==(zone), fbmc_zones)] = ref_g[t, p] / gsk4_p_plants_z[zone][t]
        else
            gsk4_matrix[t, p, findfirst(==(zone), fbmc_zones)] = 0
        end
    end
end

# GSK5: available free capacity

# denominator
gsk5_n_plants_z = Dict()
for zone in fbmc_zones
    plants_z = findall(==(zone), df_plants.zone)
    gsk5_n_plants_z[zone] = zeros(hour_count)
    for t in 1:hour_count
        num_plants_z = 0
        for p in plants_z
            if ref_g[t, p] < availability_matrix[t, p]
                num_plants_z += 1
            end
        end
        gsk5_n_plants_z[zone][t] = num_plants_z
    end
end

# nominator
gsk5_matrix = zeros(hour_count, num_p, size(fbmc_zones)[1])
for p in 1:num_p
    zone = df_plants[p, "zone"]
    for t in 1:hour_count
        if gsk5_n_plants_z[zone][t] > 0
            has_free_capacity = 0
            if ref_g[t, p] < availability_matrix[t, p]
                has_free_capacity = 1
            end
            gsk5_matrix[t, p, findfirst(==(zone), fbmc_zones)] = has_free_capacity / gsk5_n_plants_z[zone][t]
        else
            gsk5_matrix[t, p, findfirst(==(zone), fbmc_zones)] = 0
        end
    end
end

save("./flow_based_domain/gsk4_matrix_p.jld", "data", gsk4_matrix)
save("./flow_based_domain/gsk5_matrix_p.jld", "data", gsk5_matrix)