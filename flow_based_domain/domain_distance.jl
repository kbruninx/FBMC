function calculate_distance(vertices, vertices_obs)
    distance = 0
    distance_denom = 0
    for v in vertices
        distances = []
        distances_denom = []
        for v_obs in vertices_obs
            push!(distances, norm(v - v_obs))
            push!(distances_denom, norm(v_obs - [0, 0]))
        end
        distance += minimum(distances)
        distance_denom += distances_denom[findfirst(d -> d == minimum(distances), distances)]
    end
    return distance/distance_denom
end

function process_geometry(zones_s, do_central)
    println(zones_s)
    num_z = size(zones)[1]
    num_t = size(df_timestamps)[1]

    dims_to_eliminate = setdiff(1:size(zones)[1], [findfirst(z -> z == z_s, zones) for z_s in zones_s])
    eliminate_selector = ones(num_z)
    for z_s in zones_s
        eliminate_selector[findfirst(z -> z == z_s, zones)] = 0
    end

    distance_10 = []
    distance_5 = []
    distance_n = []
    #for t in 2:num_t
    for t in [400]
        println(t)
        df_ptdf_obs_t = df_ptdf_obs[df_ptdf_obs.DateTime .== df_timestamps[t], :]
        df_ptdf_10_t = df_ptdf_10[df_ptdf_10.DateTime .== df_timestamps[t], :]
        df_ptdf_5_t = df_ptdf_5[(df_ptdf_5.DateTime .== df_timestamps[t]) .& (df_ptdf_5.ram .> 0), :]
        df_ptdf_n_t = df_ptdf_n[df_ptdf_n.DateTime .== df_timestamps[t], :]

        m_obs = Model(HiGHS.Optimizer)
        m_10 = Model(HiGHS.Optimizer)
        m_5 = Model(HiGHS.Optimizer)
        m_n = Model(HiGHS.Optimizer)

        @variable(m_obs, np_obs[1:num_z])
        @variable(m_10, np_10[1:num_z])
        @variable(m_5, np_5[1:num_z])
        @variable(m_n, np_n[1:num_z])

        for j in 1:size(df_ptdf_obs_t)[1]
            @constraint(m_obs, sum(df_ptdf_obs_t[j, zones[z]]*np_obs[z] for z in 1:num_z) <= df_ptdf_obs_t[j, :ram])
        end

        for j in 1:size(df_ptdf_10_t)[1]
            @constraint(m_10, sum(df_ptdf_10_t[j, zones[z]]*np_10[z] for z in 1:num_z) <= df_ptdf_10_t[j, :ram])
        end

        for j in 1:size(df_ptdf_5_t)[1]
            @constraint(m_5, sum(df_ptdf_5_t[j, zones[z]]*np_5[z] for z in 1:num_z) <= df_ptdf_5_t[j, :ram])
        end

        for j in 1:size(df_ptdf_n_t)[1]
            @constraint(m_n, sum(df_ptdf_n_t[j, zones[z]]*np_n[z] for z in 1:num_z) <= df_ptdf_n_t[j, :ram])
        end

        poly_obs = polyhedron(m_obs, CDDLib.Library(:exact))
        poly_10 = polyhedron(m_10, CDDLib.Library(:exact))
        poly_5 = polyhedron(m_5, CDDLib.Library(:exact))
        poly_n = polyhedron(m_n, CDDLib.Library(:exact))

        poly_obs = Polyhedra.fixandeliminate(poly_obs, dims_to_eliminate, rationalize.(zeros(num_z - 2)))
        poly_10 = Polyhedra.fixandeliminate(poly_10, dims_to_eliminate, rationalize.(zeros(num_z - 2)))
        poly_5 = Polyhedra.fixandeliminate(poly_5, dims_to_eliminate, rationalize.(zeros(num_z - 2)))
        poly_n = Polyhedra.fixandeliminate(poly_n, dims_to_eliminate, rationalize.(zeros(num_z - 2)))

        plt = Plots.plot(poly_10, ratio=:equal, alpha=0.5, label=L"$\alpha = 10$")
        Plots.plot!(poly_5, ratio=:equal, alpha=0.5, label=L"$\alpha = 5$")
        Plots.plot!(poly_n, ratio=:equal, alpha=0.5, label="naive")
        Plots.plot!(poly_obs, 
            ratio=:equal, alpha=0.5, 
            label="observation", legend=true,
            title="Flow-based domain projected to the "*zones_s[1]*"-"*zones_s[2]*" plane", 
            titlefont = font(10,"Computer Modern"),
            xlabel=zones_s[2]*" net position [MW]",
            ylabel=zones_s[1]*" net position [MW]",
            xguidefontsize=9,
            yguidefontsize=9 ,
            #size=(800,600)
        )
        display(plt)
        savefig("./figures/flow_based_domain/planes/"*zones_s[1]*"_"*zones_s[2]*"_plane.png")

        try
            if do_central
                com_obs = center_of_mass(poly_obs)
                denom_com = norm(com_obs - [0,0])
                d10 = norm(center_of_mass(poly_10) - com_obs) / denom_com
                d5 = norm(center_of_mass(poly_5) - com_obs) / denom_com
                dn = norm(center_of_mass(poly_n) - com_obs) / denom_com
            else
                v_obs = points(vrep(poly_obs))
                d10 = calculate_distance(points(vrep(poly_10)), v_obs)
                d5 = calculate_distance(points(vrep(poly_5)), v_obs)
                dn = calculate_distance(points(vrep(poly_n)), v_obs)
            end

            if !isinf(d10)
                push!(distance_10, d10)
            end

            if !isinf(d5)
                push!(distance_5, d5)
            end

            if !isinf(dn)
                push!(distance_n, dn)
            end
        catch e
            println("ERROR ENCOUNTERED (not recording)")
        end
    end

    return [mean(distance_10), mean(distance_5), mean(distance_n)]
end

borders = [
    ["AT", "CZ"], ["AT", "DE_LU"], ["AT", "HU"], ["AT", "SI"], ["BE", "FR"],
    ["BE", "NL"], ["CZ", "DE_LU"], ["CZ", "PL"], ["CZ", "SK"], ["DE_LU", "FR"],
    ["DE_LU", "NL"], ["HR", "HU"], ["HR", "SI"], ["HU", "SK"], ["HU", "RO"],
    ["BE", "DE_LU"], ["DE_LU", "PL"], ["HU", "SI"], ["PL", "SK"], ["AT", "SK"],
]

border_distance_matrix = zeros(size(borders)[1], 3)
border_distance_c_matrix = zeros(size(borders)[1], 3)

for b in 1:size(borders)[1]
#for b in borders_again_index
    #border_distance_matrix[b, :] = process_geometry(borders[b], false)
    border_distance_c_matrix[b, :] = process_geometry(borders[b], true)
end

#save("./flow_based_domain/border_distance_matrix.jld", "data", border_distance_matrix)
#border_distance_matrix

#save("./flow_based_domain/border_distance_c_matrix.jld", "data", border_distance_c_matrix)
border_distance_c_matrix