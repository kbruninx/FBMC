num_p = size(df_plants)[1]
num_n = size(H_mat)[2]
num_l = size(H_mat)[1]
num_z = 12
num_t = 24 # running base case on a daily basis

function plants_at_node(n)
    node_id = df_substation_node_map[df_substation_node_map.node .== n - 1, :node_id][1]
    return findall(==(node_id), df_plants.node)
end

function nodes_in_zone(zone)
    node_ids = df_substations[df_substations.zone .== zone, :node_id]
    nodes = df_substation_node_map[in.(df_substation_node_map.node_id, [Set(node_ids)]), :node]
    return (nodes .+ 1)
end

function plants_zone_tech(zone, tech)
    return [1, 2, 3]
end

line_cap = zeros(num_l)

for l = 1:size(H_mat)[1]
    line_id = df_line_edge_map[df_line_edge_map.edge .== l - 1, :line_id][1]
    voltage = df_grid[df_grid.line_id .== line_id, :voltage][1]
    current = df_grid[df_grid.line_id .== line_id, :max_current][1]
    
    if ismissing(current)
        current = 10000
    end

    line_cap[l] = sqrt(3) * voltage * current * 2 / 1000
end

start_date = DateTime(2023, 2, 1)
end_date = DateTime(2023, 2, 2)

num_t_passed = findfirst(t->t == start_date, df_timestamps.DateTime)
day_count = Dates.value(convert(Dates.Day, end_date - start_date))

day = 1
#for day in 1:day_count
    println("DAY ", day)

    demand = zeros(num_t, num_n)

    demand_t_z_day = demand_t_z[(num_t_passed+1):(num_t_passed)+num_t, :]
    ren_gen_t_z_day = ren_gen_t_z[(num_t_passed+1):(num_t_passed)+num_t, :]
    residual_demand_t_z_day = demand_t_z_day - ren_gen_t_z_day

    for zone in zones
        z = findall(zones .== zone)[1]
        num_nodes = size(nodes_in_zone(zone))[1]
        for n in nodes_in_zone(zone) 
            demand[:, n] = residual_demand_t_z_day[:, z]/num_nodes # divide it evenly for now
        end
    end

    model = Model(HiGHS.Optimizer)
    set_optimizer_attribute(model, "presolve", "on")

    @variable(model, g[1:num_t, 1:num_p] >= 0.0)
    @variable(model, inj[1:num_t, 1:num_n])
    @variable(model, delta[1:num_t, 1:num_n])
    @variable(model, flow[1:num_t, 1:num_l])

    @constraint(model, balance_c[t=1:num_t, n=1:num_n], demand[t, n] + inj[t, n] == sum(g[t, p] for p in plants_at_node(n)))

    @constraint(model, flow_c[t=1:num_t, l=1:num_l], flow[t, l] == sum(H_mat[l, n] * delta[t, n] for n=1:num_n))
    @constraint(model, inj_c[t=1:num_t, n=1:num_n], inj[t, n] == sum(L_mat[n, q] * delta[t, q] for q=1:num_n))
    @constraint(model, cap_c_max[t=1:num_t, l=1:num_l], flow[t, l] <= line_cap[l])
    @constraint(model, cap_c_min[t=1:num_t, l=1:num_l], flow[t, l] >= -1 * line_cap[l])

    global num_t_passed += num_t  

    optimize!(model)
#end