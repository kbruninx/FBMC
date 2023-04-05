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

num_t_passed = 0
num_t = 24
day_count = 56

price_forecasts = []

#day = 1
for day in 1:day_count
    println("DAY ", day)

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

    push!(price_forecasts, JuMP.dual.(balance))
    global num_t_passed += num_t  
end

zone_list = []
for z in 3:num_z
    prices = []
    for d in 1:day_count
        for t in 1:num_t
            push!(prices, price_forecasts[d][num_t*(z-1)+t])
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