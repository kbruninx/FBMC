function plants_at_node(n)
    node_id = df_substation_node_map[df_substation_node_map.node .== n - 1, :node_id][1]
    return findall(==(node_id), df_plants.node)
end

function nodes_in_zone(zone)
    node_ids = df_substations[(df_substations.zone .== zone) .& (in.(df_substations.node_id, [Set(df_substation_node_map.node_id)])), :node_id]
    demand_keys = df_substations[(df_substations.zone .== zone) .& (in.(df_substations.node_id, [Set(df_substation_node_map.node_id)])), :demand_key]
    nodes = df_substation_node_map[in.(df_substation_node_map.node_id, [Set(node_ids)]), :node]
    return [(nodes .+ 1) demand_keys]
end

function plants_zone_tech(zone, tech)
    return [1, 2, 3]
end

coefficients_data = load(string("./base_case/coefficients.jld"))["data"]
experiment_results_alpha = coefficients_data["alpha"]
experiment_results_beta = coefficients_data["beta"]
experiment_results_gamma = coefficients_data["gamma"]
experiment_results_objective = coefficients_data["objective"]
experiment_results_timestamp = coefficients_data["timestamps"]

line_cap = zeros(num_l)

for l = 1:size(H_mat)[1]
    line_id = df_line_edge_map[df_line_edge_map.edge .== l - 1, :line_id][1]
    voltage = df_grid[df_grid.line_id .== line_id, :voltage][1]
    reactance = df_grid[df_grid.line_id .== line_id, :reactance][1]
    current = df_grid[df_grid.line_id .== line_id, :max_current][1]
    length = df_grid[df_grid.line_id .== line_id, :length][1]
    
    if ismissing(current)
        current = 3000
    end

    if length > 100
        line_cap[l] = voltage^2/(2*reactance)
    else
        line_cap[l] = sqrt(3) * voltage * current / 1000
    end
end

function optimise_base_case(num_t_passed, current_date, optimality_problem, eps_cap_pre, use_nuts)
    coeff_t_index = findlast(t -> t < current_date, experiment_results_timestamp)

    alpha = zeros(2*num_tech)
    beta = zeros(2*num_tech)
    gamma = zeros(2*num_tech)
    for z in 3:(num_z+2)
        alpha_coeffs = zeros(num_tech)
        beta_coeffs = zeros(num_tech)
        gamma_coeffs = zeros(num_tech)

        for tech in 1:num_tech
            alpha_exp = []
            beta_exp = []
            gamma_exp = []

            for e in 1:size(experiment_results_alpha[1:coeff_t_index])[1]
                push!(alpha_exp, experiment_results_alpha[e][num_tech*(z-1)+tech])
                push!(beta_exp, experiment_results_beta[e][num_tech*(z-1)+tech])
                push!(gamma_exp, experiment_results_gamma[e][num_tech*(z-1)+tech])
            end

            alpha_exp = convert(Vector{Float64}, alpha_exp)
            beta_exp = convert(Vector{Float64}, beta_exp)
            gamma_exp = convert(Vector{Float64}, gamma_exp)

            alpha_mean = mean(alpha_exp[findall(!iszero, alpha_exp)])
            beta_mean = mean(beta_exp[findall(!iszero, beta_exp)])
            gamma_mean = mean(gamma_exp[findall(!iszero, gamma_exp)])

            if !isnan(alpha_mean)
                alpha_coeffs[tech] = alpha_mean
            end
            if !isnan(beta_mean)
                beta_coeffs[tech] = beta_mean
            end
            if !isnan(gamma_mean)
                gamma_coeffs[tech] = gamma_mean
            end
        end
        
        alpha = vcat(alpha, alpha_coeffs)
        beta = vcat(beta, beta_coeffs)
        gamma = vcat(gamma, gamma_coeffs)
    end
    alpha = alpha[2*num_tech:end]
    beta = beta[2*num_tech:end]
    gamma = gamma[2*num_tech:end]

    coal_prices = coal_prices_g[(num_t_passed+1):(num_t_passed+num_t)]
    oil_prices = oil_prices_g[(num_t_passed+1):(num_t_passed+num_t)]
    gas_prices = gas_prices_g[(num_t_passed+1):(num_t_passed+num_t)]

    demand = zeros(num_t, num_n)

    demand_t_z_day = demand_t_z[(num_t_passed+1):(num_t_passed+num_t), :]
    ren_gen_t_z_day = ren_gen_t_z[(num_t_passed+1):(num_t_passed+num_t), :]
    residual_demand_t_z_day = demand_t_z_day - ren_gen_t_z_day

    for zone in zones
        z = findall(zones .== zone)[1]
        num_nodes = size(nodes_in_zone(zone))[1]
        for node_data in eachrow(nodes_in_zone(zone))
            n = node_data[1]
            demand_key = node_data[2]
            if use_nuts
                demand[:, n] = residual_demand_t_z_day[:, z] * demand_key
            else
                demand[:, n] = residual_demand_t_z_day[:, z] / num_nodes
            end
        end
    end

    #model = Model(HiGHS.Optimizer)
    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "time_limit", 180.0)
    #set_optimizer_attribute(model, "presolve", "on")

    @variable(model, g[1:num_t, 1:num_p] >= 0.0)
    @variable(model, inj[1:num_t, 1:num_n])
    @variable(model, delta[1:num_t, 1:num_n])
    @variable(model, flow[1:num_t, 1:num_l])

    if !optimality_problem
        @variable(model, eps_cap[1:num_l] >= 0.0)
    end

    @constraint(model, g_max[t=1:num_t, p=1:num_p], g[t, p] <= availability_matrix[t, p])
    @constraint(model, balance_c[t=1:num_t, n=1:num_n], demand[t, n] + inj[t, n] == sum(g[t, p] for p in plants_at_node(n)))

    @constraint(model, flow_c[t=1:num_t, l=1:num_l], flow[t, l] == sum(H_mat[l, n] * delta[t, n] for n=1:num_n))
    @constraint(model, inj_c[t=1:num_t, n=1:num_n], inj[t, n] == sum(L_mat[n, q] * delta[t, q] for q=1:num_n))

    if optimality_problem
        @constraint(model, cap_c_max[t=1:num_t, l=1:num_l], flow[t, l] <= line_cap[l] + eps_cap_pre[l])
        @constraint(model, cap_c_min[t=1:num_t, l=1:num_l], flow[t, l] >= -1 * (line_cap[l] + eps_cap_pre[l]))
    else
        @constraint(model, cap_c_max[t=1:num_t, l=1:num_l], flow[t, l] <= line_cap[l] + eps_cap[l])
        @constraint(model, cap_c_min[t=1:num_t, l=1:num_l], flow[t, l] >= -1 * (line_cap[l] + eps_cap[l]))
    end

    objective_function = 0
    for t in 1:num_t
        for z in 1:num_z
            for tech in 1:num_tech
                g_z_tech_t = sum(g[t, findall(f -> f == true, ((df_plants.zone .== zones[z]) .& (df_plants.type .== plant_types[tech])))])
                if tech == 2 || tech == 3 || tech == 5 # coal
                    objective_function += coal_prices[t]/10 * (alpha[num_tech*(z-1)+tech] + beta[num_tech*(z-1)+tech] * g_z_tech_t + gamma[num_tech*(z-1)+tech] * g_z_tech_t * g_z_tech_t)
                elseif tech == 4 # gas
                    objective_function += gas_prices[t]/10 * (alpha[num_tech*(z-1)+tech] + beta[num_tech*(z-1)+tech] * g_z_tech_t + gamma[num_tech*(z-1)+tech] * g_z_tech_t * g_z_tech_t)
                elseif tech == 6 # oil
                    objective_function += oil_prices[t]/10 * (alpha[num_tech*(z-1)+tech] + beta[num_tech*(z-1)+tech] * g_z_tech_t + gamma[num_tech*(z-1)+tech] * g_z_tech_t * g_z_tech_t)
                else
                    objective_function += alpha[num_tech*(z-1)+tech] + (beta[num_tech*(z-1)+tech] * g_z_tech_t + gamma[num_tech*(z-1)+tech] * g_z_tech_t * g_z_tech_t)
                end
            end
        end
    end

    # adding construction year consideration to objective function
    objective_construction = sum(plant_weights[p]*g[t, p] for t in 1:num_t for p in 1:num_p)
    objective_function += objective_construction

    if optimality_problem
        @objective(model, Min, objective_function)
    else
        @objective(model, Min, sum(eps_cap))
    end

    try
        optimize!(model)

        if optimality_problem
            if termination_status(model) == OPTIMAL || termination_status(model) == LOCALLY_SOLVED || termination_status(model) == ALMOST_OPTIMAL || termination_status(model) == ALMOST_LOCALLY_SOLVED
                return JuMP.value.(flow), JuMP.value.(g)
            else
                return zeros(num_t, num_l), zeros(num_t, num_p)
            end
        else
            return ceil.(JuMP.value.(eps_cap))
        end
    catch e
        println("ERROR ENCOUNTERED")
        if optimality_problem
            return zeros(num_t, num_l), zeros(num_t, num_p)
        else
            return zeros(num_l)
        end
    end
end

function run_base_case(start_date, end_date, use_nuts)
    
    num_t_passed = findfirst(t->t == start_date, df_timestamps.DateTime)
    day_count = Dates.value(convert(Dates.Day, end_date - start_date))
    
    flow_mat = []
    g_mat = []
    
    current_date = start_date
    for day in 1:day_count
        println("DAY ", day)
    
        eps_cap_pre = optimise_base_case(num_t_passed, current_date, false, missing, use_nuts)
        flow, g = optimise_base_case(num_t_passed, current_date, true, eps_cap_pre, use_nuts)
    
        if day == 1
            flow_mat = flow
            g_mat = g
        else
            flow_mat = vcat(flow_mat, flow)
            g_mat = vcat(g_mat, g)
        end
    
        num_t_passed += num_t 
        current_date = current_date + Dates.Day(1) 
    end
    
    if use_nuts
        save("./base_case/ref_flow_nuts.jld", "data", flow_mat)
        save("./base_case/ref_g_nuts.jld", "data", g_mat) 
    else
        save("./base_case/ref_flow.jld", "data", flow_mat)
        save("./base_case/ref_g.jld", "data", g_mat) 
    end
end

start_date = DateTime(2022, 9, 11)
#end_date = DateTime(2022, 9, 25)
end_date = DateTime(2023, 3, 1)

#run_base_case(start_date, end_date, true)
run_base_case(start_date, end_date, false)