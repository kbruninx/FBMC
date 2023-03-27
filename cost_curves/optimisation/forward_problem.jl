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
day_count = 1 #56

price_forecasts = []

day = 1
#for day in 1:day_count
    println("DAY ", day)

    # NON-LINEAR PART

    coal_prices = coal_prices_g[(num_t_passed+1):(num_t_passed)+num_t]
    oil_prices = oil_prices_g[(num_t_passed+1):(num_t_passed)+num_t]
    gas_prices = gas_prices_g[(num_t_passed+1):(num_t_passed)+num_t]
    eua_prices = eua_prices_g[(num_t_passed+1):(num_t_passed)+num_t]

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
                g_max_t[num_t*num_tech*(z-1)+num_t*(tech-1)+t] = g_max_g[num_tech*(z-1)+tech] 
            end
        end
    end

    model = Model(Gurobi.Optimizer)
    set_optimizer_attribute(model, "NonConvex", 2)

    @variable(model, c[1:num_z*num_tech*num_t] >= 0)
    @variable(model, g[1:num_z*num_tech*num_t] >= 0.0)

    @variable(model, np[1:num_z*num_t])

    for z in 1:num_z
        for t in 1:num_t
            for tech in 1:num_tech
                if tech == 2 || tech == 3 || tech == 5 # coal
                    @constraint(model, c[num_t*num_tech*(z-1)+num_t*(tech-1)+t] == (alpha[num_tech*(z-1)+tech] + coal_prices[t]/10 * (beta[num_tech*(z-1)+tech] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t] + gamma[num_tech*(z-1)+tech] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t]^2)))
                elseif tech == 4 # gas
                    @constraint(model, c[num_t*num_tech*(z-1)+num_t*(tech-1)+t] == (alpha[num_tech*(z-1)+tech] + gas_prices[t]/10 * (beta[num_tech*(z-1)+tech] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t] + gamma[num_tech*(z-1)+tech] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t]^2)))
                elseif tech == 6 # oil
                    @constraint(model, c[num_t*num_tech*(z-1)+num_t*(tech-1)+t] == (alpha[num_tech*(z-1)+tech] + oil_prices[t]/10 * (beta[num_tech*(z-1)+tech] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t] + gamma[num_tech*(z-1)+tech] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t]^2)))
                else
                    @constraint(model, c[num_t*num_tech*(z-1)+num_t*(tech-1)+t] == (alpha[num_tech*(z-1)+tech] + (beta[num_tech*(z-1)+tech] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t] + gamma[num_tech*(z-1)+tech] * g[num_t*num_tech*(z-1)+num_t*(tech-1)+t]^2)))
                end
            end
        end
    end

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
    @constraint(model, balance, A_balance * vcat(g, np) .== b1_balance)
    @constraint(model, cat(B_gen, B_exchange; dims=(1)) * vcat(g, np) .<= vcat(b2_gen, b2_exchange)) # combined inequality constraint

    @constraint(model, sum_z_np(np, num_t) .== 0)

    @objective(model, Min, sum(c))

    optimize!(model)

    np_obs = JuMP.value.(np)
    g_obs_orig = JuMP.value.(g)
    g_obs = JuMP.value.(g)
    c_obs = JuMP.value.(c)
    u = zeros(num_t*num_tech*num_z)
    g_resid = zeros(num_t*num_z)
    for z in 1:num_z
        for t in 1:num_t
            for tech in 1:num_tech
                if g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t] > 1
                    u[num_t*num_tech*(z-1)+num_t*(tech-1)+t] = 1
                else
                    #g_resid[num_t*(z-1)+t] = g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t]
                    #g_obs[num_t*num_tech*(z-1)+num_t*(tech-1)+t] = 0
                end
            end
        end
    end

    g_obs_delta = sum(g_obs_orig .- g_obs)

    # LINEAR PART FOR DETERMINING SHADOW PRICES
    println("LINEAR PART FOR DETERMINING SHADOW PRICES")
    model = Model(Gurobi.Optimizer)

    @variable(model, c[1:num_z*num_tech*num_t] >= 0)

    @variable(model, g[1:num_z*num_tech*num_t] >= 0)
    #@variable(model, g_delta[1:num_z*num_tech*num_t])
    #@constraint(model, g .+ g_delta .== g_obs)
    #@constraint(model, sum(g_delta) .== g_obs_delta)
    @constraint(model, g .== g_obs)

    @variable(model, np[1:num_z*num_t])
    @constraint(model, np .== np_obs)
    @constraint(model, c .== c_obs .* u)

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

    B_gen = sparse(Matrix(I, num_z*num_tech*num_t, num_z*num_tech*num_t))
    b2_gen = g_max_t

    # primal constraints
    @constraint(model, balance, A_balance * vcat(g, np) .== b1_balance)

    @objective(model, Min, c' * g)

    optimize!(model)

    global price_forecasts = vcat(price_forecasts, JuMP.dual.(balance))
    global num_t_passed += num_t  
#end

zone_list = []
for z in 3:num_z
    prices = []
    for t in 1:num_t
        push!(prices, price_forecasts[num_t*(z-1)+t])
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

print(df_price_forecast)