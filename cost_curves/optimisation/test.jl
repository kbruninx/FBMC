objective_weights = 1 .- normalize(experiment_results_objective)

df_coeffs = []
for z in 3:num_z
    alpha_coeffs = zeros(num_tech)
    beta_coeffs = zeros(num_tech)
    gamma_coeffs = zeros(num_tech)

    for tech in 1:num_tech
        alpha_exp = []
        beta_exp = []
        gamma_exp = []

        for e in 1:size(experiment_results_alpha)[1]
            push!(alpha_exp, experiment_results_alpha[e][num_tech*(z-1)+tech])
            push!(beta_exp, experiment_results_beta[e][num_tech*(z-1)+tech])
            push!(gamma_exp, experiment_results_gamma[e][num_tech*(z-1)+tech])
        end

        alpha_exp = convert(Vector{Float64}, alpha_exp)
        beta_exp = convert(Vector{Float64}, beta_exp)
        gamma_exp = convert(Vector{Float64}, gamma_exp)

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
    
    push!(df_coeffs, DataFrames.DataFrame(alpha=alpha_coeffs, beta=beta_coeffs, gamma=gamma_coeffs))
end

XLSX.writetable("cost_coefficients_alpha_fuel.xlsx",
    "AT" => df_coeffs[1],
    "BE" => df_coeffs[2],
    "CZ" => df_coeffs[3],
    "DE_LU" => df_coeffs[4],
    "FR" => df_coeffs[5],
    "HR" => df_coeffs[6],
    "HU" => df_coeffs[7],
    "NL" => df_coeffs[8],
    "PL" => df_coeffs[9],
    "RO" => df_coeffs[10],
    "SI" => df_coeffs[11],
    "SK" => df_coeffs[12],
    overwrite=true
)