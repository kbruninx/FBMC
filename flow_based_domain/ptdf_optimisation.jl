highs = optimizer_with_attributes(HiGHS.Optimizer, "Presolve" => 1) 
ipopt = optimizer_with_attributes(Ipopt.Optimizer)
alpine = optimizer_with_attributes(Alpine.Optimizer, "nlp_solver" => ipopt, "mip_solver" => highs)
minlp_solver = optimizer_with_attributes(Juniper.Optimizer, "nl_solver" => ipopt)

function optimise_direction()

    model = Model(minlp_solver)

    @variable(model, v_s_plus[1:size(df_cnecs)[1]], Bin)
    @variable(model, v_s_minus[1:size(df_cnecs)[1]], Bin)
    @variable(model, test_plus, Bin)
    @variable(model, test_minus, Bin)
    @constraint(model, v_s_plus - v_s_minus .== 1)

    M = zeros(size(O)[2], size(O)[1])
    X = zeros(size(O)[3], size(GSK1_P)[2])
    X[1, :] .= 1

    for p in 1:size(GSK1_P)[2]
        M[:, p] = O[p, :, :] * X[:, p]
    end

    objective_function = 0

    temp_line_id = 0
    s = 1
    for cnec in eachrow(df_cnecs)[179:179]
        println("Line id: ", cnec.line_id)
        println("Contingency: ", cnec.contingency)

        if ismissing(cnec.line_id)
            cnec.line_id = temp_line_id
        else
            global temp_line_id = cnec.line_id
        end
        
        for obs in eachrow(df_ptdf[(df_ptdf.line_id .== cnec.line_id) .& (df_ptdf.contingency .== cnec.contingency), :])
            edge = df_line_edge_map[df_line_edge_map.line_id .== cnec.line_id, :edge][1]
            t = findfirst(==(obs.DateTime), df_timestamps.DateTime)
            PTDF_Z = PTDF_N_C[cnec.contingency + 1] * M * GSK1_P[t, :, :]
            for zone in zones
                zone_i = findall(zones.==zone)[1]
                println(zone)
                #global objective_function += ((v_s_plus[s] - v_s_minus[s]) * obs[zone] - PTDF_Z[edge + 1, zone_i]) ^ 2
                global objective_function += ((test_plus - test_minus) * obs[zone] - PTDF_Z[edge + 1, zone_i]) ^ 2
                println(PTDF_Z[edge + 1, zone_i])
            end
        end

        global s += 1
    end

    #@NLobjective(model, Min, objective_function)
    #optimize!(model)
end

function optimise_gsk_strategy()

    model = Model(alpine)

    @variable(model, sigma[1:num_gsk_strategy, 1:12] >= 0)
    @constraint(model, sum(sigma, dims=1) .== 1)

    M = zeros(size(O)[2], size(O)[1])
    X = zeros(size(O)[3], size(GSK1_P)[2])
    X[1, :] .= 1

    for p in 1:size(GSK1_P)[2]
        M[:, p] = O[p, :, :] * X[:, p]
    end

    objective_function = 0

    temp_line_id = 0
    s = 1
    for cnec in eachrow(df_cnecs)[179:179]
        cnec_zone_i = findall(zones.==cnec.zone)[1] - 2
        println("Line id: ", cnec.line_id)
        println("Contingency: ", cnec.contingency)

        if ismissing(cnec.line_id)
            cnec.line_id = temp_line_id
        else
            temp_line_id = cnec.line_id
        end
        
        for obs in eachrow(df_ptdf[(df_ptdf.line_id .== cnec.line_id) .& (df_ptdf.contingency .== cnec.contingency), :])
            edge = df_line_edge_map[df_line_edge_map.line_id .== cnec.line_id, :edge][1]
            t = findfirst(==(obs.DateTime), df_timestamps.DateTime)

            GSK_P = [GSK1_P[t, :, :]; GSK2_P[t, :, :]; GSK3_P[t, :, :]]
            PTDF_Z = PTDF_N_C[cnec.contingency + 1] * M * sigma[:, cnec_zone_i] * GSK_P
            
            objective_function += sum(([values(obs[15:18])...] - PTDF_Z[edge + 1, zone_i]') .^ 2)
        end

        s += 1
    end

    @NLobjective(model, Min, objective_function)
    optimize!(model)
end