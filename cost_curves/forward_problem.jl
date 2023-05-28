using Dates
using JLD
using JuMP, HiGHS
using DataFrames, XLSX
using LinearAlgebra
using Statistics
using SparseArrays
using Formatting
using StatsBase

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

function has_interconnector_ex_1(z, b)
    if all_zones[z] == borders[b][1]
        return -1
    elseif all_zones[z] == borders[b][2]
            return 1
    else
        return 0
    end
end

function has_interconnector_ex_2(z, b)
    if all_zones[z] == reverse(borders[b])[1]
        return -1
    elseif all_zones[z] == reverse(borders[b])[2]
            return 1
    else
        return 0
    end
end

function get_atc_ex_1(b, t)
    return df_atc[!, borders[b][1] * "_" * borders[b][2]][t]
end

function get_atc_ex_2(b, t)
    return df_atc[!, borders[b][2] * "_" * borders[b][1]][t]
end

scenario_names = [
    #"norm_1_duality_gap_w_atc",
    "norm_1_w_atc",
    #"norm_2_duality_gap_w_atc",
    #"norm_2_w_atc",
]

months = [
    #"november", 
    "february"
]

for scenario_name in scenario_names

    coefficients_data = load(string("coefficients_", scenario_name, ".jld"))["data"]
    experiment_results_alpha = coefficients_data["alpha"]
    experiment_results_beta = coefficients_data["beta"]
    experiment_results_gamma = coefficients_data["gamma"]
    experiment_results_objective = coefficients_data["objective"]
    experiment_results_timestamp = coefficients_data["timestamps"]
    objective_weights = 1 .- normalize(experiment_results_objective)

    for month in months

        if month == "november"
            # November
            start_date = DateTime(2022, 11, 1)
            end_date = DateTime(2022, 12, 1)
        elseif month == "february"
            # February
            start_date = DateTime(2023, 2, 1)
            end_date = DateTime(2023, 3, 1)
        end

        num_t_passed = findfirst(t->t == start_date, timestamps)
        num_t = 24
        day_count = Dates.value(convert(Dates.Day, end_date - start_date))

        price_forecasts = []
        np_forecasts = []

        current_date = start_date
        for day in 1:day_count
            println("DAY ", day)

            coeff_t_index = findlast(t -> t < current_date, experiment_results_timestamp)

            alpha = zeros(2*num_tech)
            beta = zeros(2*num_tech)
            gamma = zeros(2*num_tech)
            for z in 3:(num_z+num_z_non_fbmc)
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

                    #alpha_mean = mean(alpha_exp[findall(!iszero, alpha_exp)])
                    #beta_mean = mean(beta_exp[findall(!iszero, beta_exp)])
                    #gamma_mean = mean(gamma_exp[findall(!iszero, gamma_exp)])

                    alpha_mean = mean(alpha_exp[findall(!iszero, alpha_exp)], Weights(objective_weights[findall(!iszero, alpha_exp)]))
                    beta_mean = mean(beta_exp[findall(!iszero, beta_exp)], Weights(objective_weights[findall(!iszero, beta_exp)]))
                    gamma_mean = mean(gamma_exp[findall(!iszero, gamma_exp)], Weights(objective_weights[findall(!iszero, gamma_exp)]))

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

            ch_gen_out = vec(ch_gen_out_g[(num_t_passed+1):(num_t_passed)+num_t, :])
            gb_gen_out = vec(gb_gen_out_g[(num_t_passed+1):(num_t_passed)+num_t, :])
            es_gen_out = vec(es_gen_out_g[(num_t_passed+1):(num_t_passed)+num_t, :])
            it_gen_out = vec(it_gen_out_g[(num_t_passed+1):(num_t_passed)+num_t, :])

            g_out = vcat(albe_gen_out, alde_gen_out, at_gen_out, be_gen_out, cz_gen_out, de_gen_out, fr_gen_out, hr_gen_out, hu_gen_out, nl_gen_out, pl_gen_out, ro_gen_out, si_gen_out, sk_gen_out, ch_gen_out, gb_gen_out, es_gen_out, it_gen_out) #, dk_gen_out, no_gen_out)
            g_out = convert(Vector{Float64}, g_out)

            demand_fbmc = vec(demand_g[(num_t_passed+1):(num_t_passed)+num_t, :])
            demand_non_fbmc = vec(demand_g_non_fbmc[(num_t_passed+1):(num_t_passed)+num_t, :])
            demand = vcat(demand_fbmc, demand_non_fbmc)

            ren_gen_fbmc = vec(ren_gen_g[(num_t_passed+1):(num_t_passed)+num_t, :])
            ren_gen_non_fbmc = vec(ren_gen_g_non_fbmc[(num_t_passed+1):(num_t_passed)+num_t, :])
            ren_gen = vcat(ren_gen_fbmc, ren_gen_non_fbmc)
            ren_gen = convert(Vector{Float64}, ren_gen)

            g_max_t_fbmc = zeros(num_z*num_tech*num_t)
            for z in 1:num_z
                for t in 1:num_t
                    for tech in 1:num_tech
                        g_max_t_fbmc[num_t*num_tech*(z-1)+num_t*(tech-1)+t] = g_max_g[num_tech*(z-1)+tech] 
                    end
                end
            end

            g_max_t_non_fbmc = zeros(num_z_non_fbmc*num_tech*num_t)
            for z in 1:num_z_non_fbmc
                for t in 1:num_t
                    for tech in 1:num_tech
                        g_max_t_non_fbmc[num_t*num_tech*(z-1)+num_t*(tech-1)+t] = g_max_g_non_fbmc[num_tech*(z-1)+tech] 
                    end
                end
            end

            g_max_t = vcat(g_max_t_fbmc, g_max_t_non_fbmc)

            # GENERATING BLOCK BIDS

            max_block_amount_per_tech = 30
            min_block_width = 5

            block_bids = []
            for t in 1:num_t
                block_bids_t = []
                for z in 1:(num_z+num_z_non_fbmc)
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
                        while vol_pos < g_max_t[num_t*num_tech*(z-1)+num_t*(tech-1)+t]
                            next_pos = vol_pos + block_width
                            if next_pos > g_max_t[num_t*num_tech*(z-1)+num_t*(tech-1)+t]
                                next_pos = g_max_t[num_t*num_tech*(z-1)+num_t*(tech-1)+t]
                            end
                            calc_pos = (vol_pos + next_pos)/2
                            if tech == 2 || tech == 3 || tech == 5 # coal
                                price = coal_prices[t]/10 * (alpha[num_tech*(z-1)+tech] + beta[num_tech*(z-1)+tech] * calc_pos + gamma[num_tech*(z-1)+tech] * calc_pos^2)
                            elseif tech == 4 # gas
                                price = gas_prices[t]/10 * (alpha[num_tech*(z-1)+tech] + beta[num_tech*(z-1)+tech] * calc_pos + gamma[num_tech*(z-1)+tech] * calc_pos^2)
                            elseif tech == 6 # oil
                                price = oil_prices[t]/10 * (alpha[num_tech*(z-1)+tech] + beta[num_tech*(z-1)+tech] * calc_pos + gamma[num_tech*(z-1)+tech] * calc_pos^2)
                            else
                                price = alpha[num_tech*(z-1)+tech] + (beta[num_tech*(z-1)+tech] * calc_pos + gamma[num_tech*(z-1)+tech] * calc_pos^2)
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

            total_num_bid = sum(num_bid(t, z) for z in 1:(num_z+num_z_non_fbmc) for t in 1:num_t)

            # CLEARING

            #model = Model(Gurobi.Optimizer)
            model = Model(HiGHS.Optimizer)
            set_optimizer_attribute(model, "presolve", "on")

            @variable(model, g[1:total_num_bid] >= 0.0)
            @variable(model, np[1:num_z*num_t])
            @variable(model, atc_ex_1[1:num_atc_border*num_t] >= 0.0)
            @variable(model, atc_ex_2[1:num_atc_border*num_t] >= 0.0)

            A_balance = spzeros((num_z+num_z_non_fbmc)*num_t, total_num_bid+num_z*num_t+2*num_atc_border*num_t) # contains g and np
            prev_pos = 0
            for z in 1:(num_z+num_z_non_fbmc)
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

            prev_pos += num_z*num_t

            for z in 1:(num_z+num_z_non_fbmc)
                for b in 1:num_atc_border
                    for t in 1:num_t
                        A_balance[num_t*(z-1)+t, prev_pos+num_t*(b-1)+t] = has_interconnector_ex_1(z, b) # atc exchange direction 1
                    end
                end
            end

            prev_pos += num_atc_border*num_t

            for z in 1:(num_z+num_z_non_fbmc)
                for b in 1:num_atc_border
                    for t in 1:num_t
                        A_balance[num_t*(z-1)+t, prev_pos+num_t*(b-1)+t] = has_interconnector_ex_2(z, b) # atc exchange direction 2
                    end
                end
            end

            b1_balance = demand - ren_gen

            B_gen = sparse(cat(Matrix(I, total_num_bid, total_num_bid), spzeros(total_num_bid, num_z*num_t + 2*num_atc_border*num_t); dims=(2)))
            b2_gen = vcat(collect(block_bids[t][z][1] for z in 1:(num_z+num_z_non_fbmc) for t in 1:num_t)...)

            A_exchange = spzeros(num_t, total_num_bid+num_z*num_t+2*num_atc_border*num_t)
            prev_pos = total_num_bid
            for t in 1:num_t
                for z in 1:num_z
                    A_exchange[t, prev_pos+num_t*(z-1)+t] = 1
                end
            end

            B_exchange_temp = spzeros(num_j*num_t, num_z*num_t)
            b1_exchange = spzeros(num_t)

            ram = spzeros(num_j*num_t)
            for t in 1:num_t
                df_ptdf_h = df_ptdf[df_ptdf.DateTime .== timestamps[num_t_passed+t], :]
                for j in 1:size(df_ptdf_h)[1]
                    for z in 1:num_z
                        B_exchange_temp[num_t*(j-1)+t, num_t*(z-1)+t] = df_ptdf_h[j, fbmc_zones[z]]
                        ram[num_t*(j-1)+t] = df_ptdf_h[j, :ram]
                    end
                end
            end
            
            ram = convert(Vector{Float64}, ram)
            B_exchange_ptdf = sparse(cat(spzeros(num_j*num_t, total_num_bid), B_exchange_temp, spzeros(num_j*num_t, 2*num_atc_border*num_t); dims=(2)))
            
            B_exchange_atc = spzeros(2*num_atc_border*num_t, total_num_bid + num_z*num_t + 2*num_atc_border*num_t)
            prev_pos = total_num_bid + num_z*num_t
            
            atc_1 = zeros(num_atc_border*num_t)
            for b in 1:num_atc_border
                for t in 1:num_t
                    B_exchange_atc[num_t*(b-1)+t, prev_pos+num_t*(b-1)+t] = 1 # atc exchange direction 1
                    atc_1[num_t*(b-1)+t] = get_atc_ex_1(b, t)
                end
            end

            prev_pos += num_atc_border*num_t

            atc_2 = zeros(num_atc_border*num_t)
            for b in 1:num_atc_border
                for t in 1:num_t
                    B_exchange_atc[num_t*(b-1)+t, prev_pos+num_t*(b-1)+t] = 1 # atc exchange direction 2
                    atc_2[num_t*(b-1)+t] = get_atc_ex_2(b, t)
                end
            end

            atc_1 = convert(Vector{Float64}, atc_1)
            atc_2 = convert(Vector{Float64}, atc_2)

            b2_exchange = vcat(ram, remove_missing(atc_1), remove_missing(atc_2))
            B_exchange = sparse(cat(B_exchange_ptdf, B_exchange_atc; dims=(1)))

            # primal constraints
            @constraint(model, balance, A_balance * vcat(g, np, atc_ex_1, atc_ex_2) .== b1_balance)
            @constraint(model, cat(B_gen, B_exchange; dims=(1)) * vcat(g, np, atc_ex_1, atc_ex_2) .<= vcat(b2_gen, b2_exchange)) # combined inequality constraint

            @constraint(model, sum_z_np(np, num_t) .== 0)

            c = vcat(collect(block_bids[t][z][2] for z in 1:(num_z+num_z_non_fbmc) for t in 1:num_t)...)
            @objective(model, Min, c' * g)

            try
                optimize!(model)
                push!(price_forecasts, JuMP.dual.(balance))
                push!(np_forecasts, JuMP.value.(np))
            catch e
                push!(price_forecasts, zeros((num_z+num_z_non_fbmc)*num_t))
                push!(np_forecasts, zeros(num_z*num_t))
                println("ERROR ENCOUNTERED")
            end

            num_t_passed += num_t  
            current_date = current_date + Dates.Day(1)
        end

        zone_list = []
        np_zone_list = []

        for z in 3:(num_z+num_z_non_fbmc)
            prices = []
            for d in 1:day_count
                for t in 1:num_t
                    push!(prices, price_forecasts[d][num_t*(z-1)+t])
                end
            end
            push!(zone_list, prices)
        end

        for z in 3:num_z
            nps = []
            for d in 1:day_count
                for t in 1:num_t
                    push!(nps, np_forecasts[d][num_t*(z-1)+t])
                end
            end
            push!(np_zone_list, nps)
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
            SK=zone_list[12],
            CH=zone_list[13],
            GB=zone_list[14],
            ES=zone_list[15],
            IT_NORD=zone_list[16],
        )

        df_np_forecast = DataFrames.DataFrame(
            AT=np_zone_list[1],
            BE=np_zone_list[2],
            CZ=np_zone_list[3],
            DE_LU=np_zone_list[4],
            FR=np_zone_list[5],
            HR=np_zone_list[6],
            HU=np_zone_list[7],
            NL=np_zone_list[8],
            PL=np_zone_list[9],
            RO=np_zone_list[10],
            SI=np_zone_list[11],
            SK=np_zone_list[12],
        )

        XLSX.writetable(string("price_forecast_full_tso_", scenario_name, "_", month, ".xlsx"), df_price_forecast)
        XLSX.writetable(string("np_forecast_full_tso_", scenario_name, "_", month, ".xlsx"), df_np_forecast)

    end
end