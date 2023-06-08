using JuMP, Ipopt, Juniper

# initializing MINLP solver
ipopt = optimizer_with_attributes(Ipopt.Optimizer)
minlp_solver = optimizer_with_attributes(Juniper.Optimizer, "nl_solver" => ipopt)
model = Model(minlp_solver)

# variable for GSK selection per zone
@variable(model, sigma[1:num_gsk_strategy, 1:12] >= 0)

# relative selection preference summing up to 1
@constraint(model, sum(sigma, dims=1) .== 1)

# simple plant to node map
M = zeros(size(O)[2], size(O)[1])
X = zeros(size(O)[3], size(GSK1_P)[2])
X[1, :] .= 1

for p in 1:size(GSK1_P)[2]
    M[:, p] = O[p, :, :] * X[:, p]
end

objective_function = 0

s = 1
# gathering CNEC observations for fitting
for cnec in eachrow(df_cnecs)
    cnec_zone_i = findall(zones .== cnec.zone)[1] - 2
    
    for obs in eachrow(df_ptdf[(df_ptdf.line_id .== cnec.line_id) .& (df_ptdf.contingency .== cnec.contingency) .& (df_ptdf.DateTime .>= start_date) .& (df_ptdf.DateTime .< end_date), :])
        edge = df_line_edge_map[df_line_edge_map.line_id .== cnec.line_id, :edge][1]
        t = findfirst(==(obs.DateTime), df_timestamps[(df_timestamps.DateTime .>= start_date) .& (df_timestamps.DateTime .< end_date), :].DateTime)

        # calculating zonal PTDFs for all GSK strategies
        PTDF_Z = [
            PTDF_N_C[cnec.contingency + 1] * M * GSK1_P[t, :, :];;;
            PTDF_N_C[cnec.contingency + 1] * M * GSK2_P[t, :, :];;;
            PTDF_N_C[cnec.contingency + 1] * M * GSK3_P[t, :, :];;;
            PTDF_N_C[cnec.contingency + 1] * M * GSK4_P[t, :, :];;;
            PTDF_N_C[cnec.contingency + 1] * M * GSK5_P[t, :, :];;;
        ]

        objective_function += sum((
            [values(obs[5:18])...] - [sigma[:, cnec_zone_i]' * PTDF_Z[edge + 1, i, :] for i = 1:14]
        ) .^ 2) # L2-norm
    end

    s += 1
end

@NLobjective(model, Min, objective_function)
optimize!(model)

JuMP.value.(sigma)