using JuMP, BilevelJuMP, Gurobi, Dualization
using DataFrames, XLSX
using LinearAlgebra
using Alpine
using Ipopt
using Statistics
using QuadraticToBinary
using Plots
using SparseArrays
using Formatting
using Juniper

function sum_z_np(np, num_t)
    result = []
    for t in 1:num_t
        sum = 0
        for z in 1:num_z
            sum += np[num_t*(z-1)+t]
        end
        push!(result, sum)
    end
    return result
end

current_step = 1466 # starting after 10 days

price_forecasts = []

iteration_step = 1
while (current_step + moving_step) <= data_length
    println("DAY ", iteration_step, ": STARTING INVERSE OPTIMISATION...")
   
    # INVERSE OPTIMISATION
    num_t = training_window
    num_t_passed = current_step - training_window

    coal_prices = coal_prices_g[(num_t_passed+1):(num_t_passed)+num_t]
    oil_prices = oil_prices_g[(num_t_passed+1):(num_t_passed)+num_t]
    gas_prices = gas_prices_g[(num_t_passed+1):(num_t_passed)+num_t]
    eua_prices = eua_prices_g[(num_t_passed+1):(num_t_passed)+num_t]

    at_obs = vec(at_obs_g[(num_t_passed+1):(num_t_passed)+num_t, :])
    be_obs = vec(be_obs_g[(num_t_passed+1):(num_t_passed)+num_t, :])
    cz_obs = vec(cz_obs_g[(num_t_passed+1):(num_t_passed)+num_t, :])
    de_obs = vec(de_obs_g[(num_t_passed+1):(num_t_passed)+num_t, :])
    fr_obs = vec(fr_obs_g[(num_t_passed+1):(num_t_passed)+num_t, :])
    hr_obs = vec(hr_obs_g[(num_t_passed+1):(num_t_passed)+num_t, :])
    hu_obs = vec(hu_obs_g[(num_t_passed+1):(num_t_passed)+num_t, :])
    nl_obs = vec(nl_obs_g[(num_t_passed+1):(num_t_passed)+num_t, :])
    pl_obs = vec(pl_obs_g[(num_t_passed+1):(num_t_passed)+num_t, :])
    ro_obs = vec(ro_obs_g[(num_t_passed+1):(num_t_passed)+num_t, :])
    si_obs = vec(si_obs_g[(num_t_passed+1):(num_t_passed)+num_t, :])
    sk_obs = vec(sk_obs_g[(num_t_passed+1):(num_t_passed)+num_t, :])
    albe_obs = zeros(num_t*num_tech)
    alde_obs = zeros(num_t*num_tech)

    g_obs = vcat(albe_obs, alde_obs, at_obs, be_obs, cz_obs, de_obs, fr_obs, hr_obs, hu_obs, nl_obs, pl_obs, ro_obs, si_obs, sk_obs)
    g_obs = convert(Vector{Float64}, g_obs)

    demand = vec(demand_g[(num_t_passed+1):(num_t_passed)+num_t, :])
    lambda_obs = vec(lambda_obs_g[(num_t_passed+1):(num_t_passed)+num_t, :])
    np_obs = vec(np_obs_g[(num_t_passed+1):(num_t_passed)+num_t, :])
    demand = convert(Vector{Float64}, demand)
    lambda_obs = convert(Vector{Float64}, lambda_obs)
    np_obs = convert(Vector{Float64}, np_obs)

    ptdf_z = vec(ptdf_z_g[(num_t_passed*num_j+1):(num_t_passed*num_j)+num_t*num_j, :])
    ram = vec(ram_g[(num_t_passed*num_j+1):(num_t_passed*num_j)+num_t*num_j, :])
    ptdf_z = convert(Vector{Float64}, ptdf_z)
    ram = convert(Vector{Float64}, ram)

    ren_gen = vec(ren_gen_g[(num_t_passed+1):(num_t_passed)+num_t, :])
    ren_gen = convert(Vector{Float64}, ren_gen)

    g_max_t = zeros(num_z*num_tech*num_t)
    for z in 1:num_z
        for t in 1:num_t
            for tech in 1:num_tech
                g_max_t[num_t*num_tech*(z-1)+num_t*(tech-1)+t] = g_max_g[num_tech*(z-1)+tech] 
            end
        end
    end

    model = direct_model(Gurobi.Optimizer())
    set_time_limit_sec(model, 120.0)
    set_optimizer_attribute(model, "NonConvex", 2)

    @variable(model, c[1:num_z*num_tech*num_t] >= 0)

    @variable(model, alpha[1:num_z*num_tech] >= 0)
    @variable(model, beta[1:num_z*num_tech] >= 0)
    @variable(model, gamma[1:num_z*num_tech] >= 0)

    for z in 1:num_z
        for tech in 3:num_tech
            if tech == 7
                @constraint(model, gamma[num_tech*(z-1)+tech] == 0)
            end
        end
        for t in 1:num_t
            for tech in 1:num_tech

                if tech == 1 || tech == 8 || tech == 9 || tech == 10 # without fuel-responsiveness
                    @constraint(model, c[num_t*num_tech*(z-1)+num_t*(tech-1)+t] == alpha[num_tech*(z-1)+tech] + (beta[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t] + gamma[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t]^2))
                elseif tech == 2 || tech == 3 || tech == 5 # coal
                    @constraint(model, c[num_t*num_tech*(z-1)+num_t*(tech-1)+t] == coal_prices[t]/10 * (alpha[num_tech*(z-1)+tech] + beta[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t] + gamma[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t]^2))
                elseif tech == 4 # gas
                    @constraint(model, c[num_t*num_tech*(z-1)+num_t*(tech-1)+t] == gas_prices[t]/10 * (alpha[num_tech*(z-1)+tech] + beta[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t] + gamma[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t]^2))
                elseif tech == 6 # oil
                    @constraint(model, c[num_t*num_tech*(z-1)+num_t*(tech-1)+t] == oil_prices[t]/10 * (alpha[num_tech*(z-1)+tech] + beta[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t] + gamma[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t]^2))
                elseif tech == 7 # hydro
                    @constraint(model, c[num_t*num_tech*(z-1)+num_t*(tech-1)+t] == alpha[num_tech*(z-1)+tech] + beta[num_tech*(z-1)+tech] * g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t])
                end
            end
        end
    end

    @variable(model, eps1[1:num_z*num_t] >= 0)
    @variable(model, eps2[1:num_z*num_t] >= 0)
    #@variable(model, eps3[1:num_z*num_t] >= 0)
    #@variable(model, eps4[1:num_z*num_t] >= 0)

    @variable(model, g[1:num_z*num_tech*num_t] >= 0)
    @variable(model, np[1:num_z*num_t])

    # dual variables
    @variable(model, lambda[1:num_z*num_t])
    @variable(model, lambda_exchange[1:num_t])

    @variable(model, mu_gen[1:num_z*num_tech*num_t] <= 0)
    @variable(model, mu_exchange[1:num_j*num_t] <= 0)

    @constraint(model, g .== g_obs) # observe generation

    A_balance = spzeros(num_z*num_t, num_z*num_tech*num_t+num_z*num_t) # contains g and np
    prev_pos = num_z*num_tech*num_t
    for z in 1:num_z
        for t in 1:num_t
            for tech in 1:num_tech
                A_balance[num_t*(z-1)+t, num_t*num_tech*(z-1)+num_t*(tech-1)+t] = 1
            end
            A_balance[num_t*(z-1)+t, prev_pos+num_t*(z-1)+t] = -1 # np
        end
    end

    b1_balance = demand - ren_gen

    B_gen = sparse(cat(Matrix(I, num_z*num_tech*num_t, num_z*num_tech*num_t), spzeros(num_z*num_tech*num_t, num_z*num_t); dims=(2)))
    b2_gen = g_max_t

    A_exchange = spzeros(num_t, num_z*num_tech*num_t+num_z*num_t)
    for t in 1:num_t
        for z in 1:num_z
            A_exchange[t, prev_pos+num_t*(z-1)+t] = 1
        end
    end

    b1_exchange = spzeros(num_t)

    B_exchange_temp = spzeros(num_j*num_t, num_z*num_t)
    b2_exchange = ram
    for t in 1:num_t
        for j in 1:num_j
            for z in 1:num_z
                B_exchange_temp[num_t*(j-1)+t, num_t*(z-1)+t] = ptdf_z[num_j*num_t*(z-1) + num_j*(t-1) + j]
            end
        end
    end
    B_exchange = sparse(cat(spzeros(num_j*num_t, num_z*num_tech*num_t), B_exchange_temp; dims=(2)))

    # primal constraints
    #@constraint(model, balance, A_balance * vcat(g, np) .== b1_balance) # demand equality
    @constraint(model, B_exchange * vcat(g, np) .<= b2_exchange) # net position, rams and ptdfs
    @constraint(model, sum_z_np(np, num_t) .== 0)

    # dual constraints
    @constraint(model, cat(A_balance, A_exchange; dims=(1))' * vcat(lambda, lambda_exchange) .+ cat(B_gen, B_exchange; dims=(1))' * vcat(mu_gen, mu_exchange) .== vcat(c, spzeros(num_z*num_t)))

    # strong duality gap theorem
    @constraint(model, sum(c) .- b1_balance' * lambda .- b1_exchange' * lambda_exchange .- b2_gen' * mu_gen .- b2_exchange' * mu_exchange == 0)

    @constraint(model, lambda .- lambda_obs .== eps1 .- eps2)
    #@constraint(model, np .- np_obs .== eps3 .- eps4)

    u = ones(num_z*num_t)
    @objective(model, Min, eps1' * u + eps2' * u)
    #@objective(model, Min, (eps1' * u + eps2' * u)/(num_z*num_t*maximum(lambda_obs)) + (eps3' * u + eps4' * u)/(num_z*num_t*maximum(np_obs)))

    optimize!(model)

    if termination_status(model) == OPTIMAL 
        println("DAY ", iteration_step, ": OPTIMAL IO SOLUTION FOUND")
        alpha_r = JuMP.value.(alpha)
        beta_r = JuMP.value.(beta)
        gamma_r = JuMP.value.(gamma)

        # MERGING COST COEFFICIENT RESULTS
        

        println("DAY ", iteration_step, ": STARTING FORWARD PROBLEM...")
        # FORWARD PROBLEM
        num_t_passed = current_step
        num_t = moving_step

        coal_prices = coal_prices_g[(num_t_passed+1):(num_t_passed)+num_t]
        oil_prices = oil_prices_g[(num_t_passed+1):(num_t_passed)+num_t]
        gas_prices = gas_prices_g[(num_t_passed+1):(num_t_passed)+num_t]
        eua_prices = eua_prices_g[(num_t_passed+1):(num_t_passed)+num_t]

        at_gen_out = vec(at_gen_out_g[(num_t_passed+1):(num_t_passed)+num_t, :])
        be_gen_out = vec(be_gen_out_g[(num_t_passed+1):(num_t_passed)+num_t, :])
        cz_gen_out = vec(cz_gen_out_g[(num_t_passed+1):(num_t_passed)+num_t, :])
        de_gen_out = vec(de_gen_out_g[(num_t_passed+1):(num_t_passed)+num_t, :])
        fr_gen_out = vec(fr_gen_out_g[(num_t_passed+1):(num_t_passed)+num_t, :])
        hr_gen_out = vec(hr_gen_out_g[(num_t_passed+1):(num_t_passed)+num_t, :])
        hu_gen_out = vec(hu_gen_out_g[(num_t_passed+1):(num_t_passed)+num_t, :])
        nl_gen_out = vec(nl_gen_out_g[(num_t_passed+1):(num_t_passed)+num_t, :])
        pl_gen_out = vec(pl_gen_out_g[(num_t_passed+1):(num_t_passed)+num_t, :])
        ro_gen_out = vec(ro_gen_out_g[(num_t_passed+1):(num_t_passed)+num_t, :])
        si_gen_out = vec(si_gen_out_g[(num_t_passed+1):(num_t_passed)+num_t, :])
        sk_gen_out = vec(sk_gen_out_g[(num_t_passed+1):(num_t_passed)+num_t, :])
        albe_gen_out = zeros(num_t*num_tech)
        alde_gen_out = zeros(num_t*num_tech)

        g_out = vcat(albe_gen_out, alde_gen_out, at_gen_out, be_gen_out, cz_gen_out, de_gen_out, fr_gen_out, hr_gen_out, hu_gen_out, nl_gen_out, pl_gen_out, ro_gen_out, si_gen_out, sk_gen_out)
        g_out = convert(Vector{Float64}, g_out)

        demand = vec(demand_g[(num_t_passed+1):(num_t_passed)+num_t, :])
        demand = convert(Vector{Float64}, demand)

        ptdf_z = vec(ptdf_z_g[(num_t_passed*num_j+1):(num_t_passed*num_j)+num_t*num_j, :])
        ram = vec(ram_g[(num_t_passed*num_j+1):(num_t_passed*num_j)+num_t*num_j, :])
        ptdf_z = convert(Vector{Float64}, ptdf_z)
        ram = convert(Vector{Float64}, ram)

        ren_gen = vec(ren_gen_g[(num_t_passed+1):(num_t_passed)+num_t, :])
        ren_gen = convert(Vector{Float64}, ren_gen)

        g_max_t = zeros(num_z*num_tech*num_t)
        for z in 1:num_z
            for t in 1:num_t
                for tech in 1:num_tech
                    g_max_t[num_t*num_tech*(z-1)+num_t*(tech-1)+t] = g_max_g[num_tech*(z-1)+tech] - g_out[num_t*num_tech*(z-1)+num_t*(tech-1)+t] 
                end
            end
        end

        # GENERATING BLOCK BIDS

        max_block_amount_per_tech = 30
        min_block_width = 5

        block_bids = []
        for t in 1:num_t
            block_bids_t = []
            for z in 1:num_z
                block_bids_zone_volumes = []
                block_bids_zone_prices = []
                for tech in 1:num_tech
                    prices = []
                    volumes = []
                    positions = []

                    block_width = floor(g_max_t[num_t*num_tech*(z-1)+num_t*(tech-1)+t]/max_block_amount_per_tech)
                    if block_width < 5
                        block_width = 5
                    end

                    vol_pos = 0
                    while vol_pos < g_max_g[num_tech*(z-1)+tech]
                        next_pos = vol_pos + block_width
                        if next_pos > g_max_g[num_tech*(z-1)+tech]
                            next_pos = g_max_g[num_tech*(z-1)+tech]
                        end
                        calc_pos = (vol_pos + next_pos)/2
                        if tech == 2 || tech == 3 || tech == 5 # coal
                            price = coal_prices[t]/10 * (alpha_r[num_tech*(z-1)+tech] + beta_r[num_tech*(z-1)+tech] * calc_pos + gamma_r[num_tech*(z-1)+tech] * calc_pos^2)
                        elseif tech == 4 # gas
                            price = gas_prices[t]/10 * (alpha_r[num_tech*(z-1)+tech] + beta_r[num_tech*(z-1)+tech] * calc_pos + gamma_r[num_tech*(z-1)+tech] * calc_pos^2)
                        elseif tech == 6 # oil
                            price = oil_prices[t]/10 * (alpha_r[num_tech*(z-1)+tech] + beta_r[num_tech*(z-1)+tech] * calc_pos + gamma_r[num_tech*(z-1)+tech] * calc_pos^2)
                        else
                            price = alpha_r[num_tech*(z-1)+tech] + (beta_r[num_tech*(z-1)+tech] * calc_pos + gamma_r[num_tech*(z-1)+tech] * calc_pos^2)
                        end
                        push!(volumes, next_pos - vol_pos)
                        push!(prices, price)
                        push!(positions, calc_pos)
                        vol_pos = next_pos
                    end
                    block_bids_zone_volumes = vcat(block_bids_zone_volumes, volumes)
                    block_bids_zone_prices = vcat(block_bids_zone_prices, prices)
                end
                push!(block_bids_t, [block_bids_zone_volumes, block_bids_zone_prices])
            end
            push!(block_bids, block_bids_t)
        end

        function num_bid(t, z)
            return size(block_bids[t][z][1])[1]
        end

        total_num_bid = sum(num_bid(t, z) for z in 1:num_z for t in 1:num_t)

        # CLEARING

        model = Model(Gurobi.Optimizer)

        @variable(model, g[1:total_num_bid] >= 0.0)
        @variable(model, np[1:num_z*num_t])

        A_balance = spzeros(num_z*num_t, total_num_bid+num_z*num_t) # contains g and np
        prev_pos = 0
        for z in 1:num_z
            for t in 1:num_t
                for bid in 1:num_bid(t, z)
                    A_balance[num_t*(z-1)+t, prev_pos + 1] = 1
                    prev_pos += 1
                end
            end
        end

        for z in 1:num_z
            for t in 1:num_t
                A_balance[num_t*(z-1)+t, prev_pos+num_t*(z-1)+t] = -1 # np
            end
        end

        b1_balance = demand - ren_gen

        B_gen = sparse(cat(Matrix(I, total_num_bid, total_num_bid), spzeros(total_num_bid, num_z*num_t); dims=(2)))
        b2_gen = vcat(collect(block_bids[t][z][1] for z in 1:num_z for t in 1:num_t)...)

        A_exchange = spzeros(num_t, total_num_bid+num_z*num_t)
        for t in 1:num_t
            for z in 1:num_z
                A_exchange[t, prev_pos+num_t*(z-1)+t] = 1
            end
        end

        b1_exchange = spzeros(num_t)

        B_exchange_temp = spzeros(num_j*num_t, num_z*num_t)
        b2_exchange = ram
        for t in 1:num_t
            for j in 1:num_j
                for z in 1:num_z
                    B_exchange_temp[num_t*(j-1)+t, num_t*(z-1)+t] = ptdf_z[num_j*num_t*(z-1) + num_j*(t-1) + j]
                end
            end
        end
        B_exchange = sparse(cat(spzeros(num_j*num_t, total_num_bid), B_exchange_temp; dims=(2)))

        # primal constraints
        @constraint(model, balance, A_balance * vcat(g, np) .== b1_balance)
        @constraint(model, cat(B_gen, B_exchange; dims=(1)) * vcat(g, np) .<= vcat(b2_gen, b2_exchange)) # combined inequality constraint

        @constraint(model, sum_z_np(np, num_t) .== 0)

        c = vcat(collect(block_bids[t][z][2] for z in 1:num_z for t in 1:num_t)...)
        @objective(model, Min, c' * g)
        optimize!(model)

        println("DAY ", iteration_step, ": FORECAST DONE.")
        """
        compute_conflict!(model)
        list_of_conflicting_constraints = ConstraintRef[]
        for (F, S) in list_of_constraint_types(model)
            for con in all_constraints(model, F, S)
                if MOI.get(model, MOI.ConstraintConflictStatus(), con) == MOI.IN_CONFLICT
                    push!(list_of_conflicting_constraints, con)
                end
            end
        end
        println(list_of_conflicting_constraints)
        """

        push!(price_forecasts, JuMP.dual.(balance))
    else
        println("DAY ", iteration_step, ": NO IO SOLUTION FOUND, SKIPPING.")
        push!(price_forecasts, zeros(num_z*num_t))
    end

    global current_step += moving_step
    global iteration_step += 1
end

zone_list = []
for z in 3:num_z
    prices = []
    for d in 1:trunc(Int, (data_length-training_window)/24)
        for t in 1:moving_step
            push!(prices, price_forecasts[d][moving_step*(z-1)+t])
        end
    end
    push!(zone_list, prices)
end

df_price_forecast = DataFrames.DataFrame(
    AT=zone_list[1],
    BE=zone_list[2],
    CZ=zone_list[3],
    DE_LU=zone_list[4],
    FR=zone_list[5],
    HR=zone_list[6],
    HU=zone_list[7],
    NL=zone_list[8],
    PL=zone_list[9],
    RO=zone_list[10],
    SI=zone_list[11],
    SK=zone_list[12]
)