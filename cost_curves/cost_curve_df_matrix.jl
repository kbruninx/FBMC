
scenario_names = [
    "norm_1_duality_gap_w_atc",
    "norm_1_w_atc",
    "norm_2_duality_gap_w_atc",
    "norm_2_w_atc",
]

alpha_scenario = zeros(4,(num_z+num_z_non_fbmc),num_tech)
beta_scenario = zeros(4,(num_z+num_z_non_fbmc),num_tech)
gamma_scenario = zeros(4,(num_z+num_z_non_fbmc),num_tech)

for s in 1:size(scenario_names)[1]
    coefficients_data = load(string("./cost_curves/coefficients_", scenario_names[s], ".jld"))["data"]
    experiment_results_alpha = coefficients_data["alpha"]
    experiment_results_beta = coefficients_data["beta"]
    experiment_results_gamma = coefficients_data["gamma"]
    experiment_results_objective = coefficients_data["objective"]
    objective_weights = 1 .- normalize(experiment_results_objective)

    num_z = 14 # including ALBE and ALDE
    num_z_non_fbmc = 4
    num_tech = 10

    alpha = zeros((num_z+num_z_non_fbmc),num_tech)
    beta = zeros((num_z+num_z_non_fbmc),num_tech)
    gamma = zeros((num_z+num_z_non_fbmc),num_tech)

    for z in 3:(num_z+num_z_non_fbmc)
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
        
        alpha[z, :] = alpha_coeffs
        beta[z, :] = beta_coeffs
        gamma[z, :] = gamma_coeffs
    end

    alpha_scenario[s, :, :] = alpha
    beta_scenario[s, :, :] = beta
    gamma_scenario[s, :, :] = gamma
end

save("./cost_curves/coeffs_matrix.jld", "alpha", alpha_scenario, "beta", beta_scenario, "gamma", gamma_scenario)
