sigma = load("./flow_based_domain/sigma.jld")["data"]
sigma[sigma .< 0.01] .= 0

for z in 1:12
    sigma[:, z] = sigma[:, z] ./ sum(sigma[:, z])
end

save("./flow_based_domain/sigma.jld", "data", sigma)