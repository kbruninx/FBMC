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

function remove_missing(mx)
    mx[ismissing.(mx)] .= 0.0
    return mx
end

zones =  ["ALBE", "ALDE", "AT", "BE", "CZ", "DE_LU", "FR", "HR", "HU", "NL", "PL", "RO", "SI", "SK"]

plant_types = [
    "biomass", "brown_coal", "coal_gas", "natural_gas", "hard_coal", "oil", "hydro", 
    "nuclear", "waste", "other"
] 

xf_atc = XLSX.readxlsx("./data/validation/atc.xlsx")
xf_da_prices = XLSX.readxlsx("./data/validation/day_ahead_prices.xlsx")
xf_da_prices_non_fbmc = XLSX.readxlsx("./data/validation/day_ahead_prices_non_fbm.xlsx")
xf_demand = XLSX.readxlsx("./data/validation/demand.xlsx")
xf_fuel_prices = XLSX.readxlsx("./data/validation/fuel_prices.xlsx")
xf_generation = XLSX.readxlsx("./data/validation/generation.xlsx")
xf_ren_gen = XLSX.readxlsx("./data/validation/renewable_generation.xlsx")
xf_netpos = XLSX.readxlsx("./data/validation/net_positions.xlsx")
xf_ptdf = XLSX.readxlsx("./data/validation/ptdfs.xlsx")
xf_capacities = XLSX.readxlsx("./data/installed_capacities_corrected.xlsx")

num_z = 14 # including ALBE and ALDE
num_val_t = 1346 # size of validation set
num_tech = 10
num_j = 133 # maximum amont on CNEs at a given time in the validation dataset

coal_prices_g = vec(remove_missing(xf_fuel_prices["Sheet1"][sprintf1("B2:B%d", num_val_t+1)]))
oil_prices_g = vec(remove_missing(xf_fuel_prices["Sheet1"][sprintf1("C2:C%d", num_val_t+1)]))
gas_prices_g = vec(remove_missing(xf_fuel_prices["Sheet1"][sprintf1("D2:D%d", num_val_t+1)]))
eua_prices_g = vec(remove_missing(xf_fuel_prices["Sheet1"][sprintf1("E2:E%d", num_val_t+1)]))

coal_prices_g = convert(Vector{Float64}, coal_prices_g)
oil_prices_g = convert(Vector{Float64}, oil_prices_g)
gas_prices_g = convert(Vector{Float64}, gas_prices_g)
eua_prices_g = convert(Vector{Float64}, eua_prices_g)

at_obs_g = remove_missing(xf_generation["AT"][sprintf1("B2:K%d", num_val_t+1)])
be_obs_g = remove_missing(xf_generation["BE"][sprintf1("B2:K%d", num_val_t+1)])
cz_obs_g = remove_missing(xf_generation["CZ"][sprintf1("B2:K%d", num_val_t+1)])
de_obs_g = remove_missing(xf_generation["DE_LU"][sprintf1("B2:K%d", num_val_t+1)])
fr_obs_g = remove_missing(xf_generation["FR"][sprintf1("B2:K%d", num_val_t+1)])
hr_obs_g = remove_missing(xf_generation["HR"][sprintf1("B2:K%d", num_val_t+1)])
hu_obs_g = remove_missing(xf_generation["HU"][sprintf1("B2:K%d", num_val_t+1)])
nl_obs_g = remove_missing(xf_generation["NL"][sprintf1("B2:K%d", num_val_t+1)])
pl_obs_g = remove_missing(xf_generation["PL"][sprintf1("B2:K%d", num_val_t+1)])
ro_obs_g = remove_missing(xf_generation["RO"][sprintf1("B2:K%d", num_val_t+1)])
si_obs_g = remove_missing(xf_generation["SI"][sprintf1("B2:K%d", num_val_t+1)])
sk_obs_g = remove_missing(xf_generation["SK"][sprintf1("B2:K%d", num_val_t+1)])

demand_g = remove_missing(xf_demand["Sheet1"][sprintf1("B2:O%d", num_val_t+1)]) # [z+t]
lambda_obs_g = remove_missing(xf_da_prices["Sheet1"][sprintf1("B2:O%d", num_val_t+1)]) # [z+t]
np_obs_g = remove_missing(xf_netpos["Sheet1"][sprintf1("B2:O%d", num_val_t+1)]) # [z+t]

ptdf_z_g = remove_missing(xf_ptdf["Sheet1"][sprintf1("D2:Q%d", num_val_t*num_j+1)])
ram_g = remove_missing(xf_ptdf["Sheet1"][sprintf1("C2:C%d", num_val_t*num_j+1)])

g_max_g = vec(xf_capacities["Sheet1"]["B2:O11"]) # [z+tech]

ren_gen_g = remove_missing(xf_ren_gen["Sheet1"][sprintf1("B2:O%d", num_val_t+1)]) # [z+t]

xf_cost_coeffs = XLSX.readxlsx("./cost_coefficients_with_np.xlsx")

at_a = xf_cost_coeffs["AT"]["A2:A11"]
be_a = xf_cost_coeffs["BE"]["A2:A11"]
cz_a = xf_cost_coeffs["CZ"]["A2:A11"]
de_a = xf_cost_coeffs["DE_LU"]["A2:A11"]
fr_a = xf_cost_coeffs["FR"]["A2:A11"]
hr_a = xf_cost_coeffs["HR"]["A2:A11"]
hu_a = xf_cost_coeffs["HU"]["A2:A11"]
nl_a = xf_cost_coeffs["NL"]["A2:A11"]
pl_a = xf_cost_coeffs["PL"]["A2:A11"]
ro_a = xf_cost_coeffs["RO"]["A2:A11"]
si_a = xf_cost_coeffs["SI"]["A2:A11"]
sk_a = xf_cost_coeffs["SK"]["A2:A11"]

at_b = xf_cost_coeffs["AT"]["B2:B11"]
be_b = xf_cost_coeffs["BE"]["B2:B11"]
cz_b = xf_cost_coeffs["CZ"]["B2:B11"]
de_b = xf_cost_coeffs["DE_LU"]["B2:B11"]
fr_b = xf_cost_coeffs["FR"]["B2:B11"]
hr_b = xf_cost_coeffs["HR"]["B2:B11"]
hu_b = xf_cost_coeffs["HU"]["B2:B11"]
nl_b = xf_cost_coeffs["NL"]["B2:B11"]
pl_b = xf_cost_coeffs["PL"]["B2:B11"]
ro_b = xf_cost_coeffs["RO"]["B2:B11"]
si_b = xf_cost_coeffs["SI"]["B2:B11"]
sk_b = xf_cost_coeffs["SK"]["B2:B11"]

at_g = xf_cost_coeffs["AT"]["C2:C11"]
be_g = xf_cost_coeffs["BE"]["C2:C11"]
cz_g = xf_cost_coeffs["CZ"]["C2:C11"]
de_g = xf_cost_coeffs["DE_LU"]["C2:C11"]
fr_g = xf_cost_coeffs["FR"]["C2:C11"]
hr_g = xf_cost_coeffs["HR"]["C2:C11"]
hu_g = xf_cost_coeffs["HU"]["C2:C11"]
nl_g = xf_cost_coeffs["NL"]["C2:C11"]
pl_g = xf_cost_coeffs["PL"]["C2:C11"]
ro_g = xf_cost_coeffs["RO"]["C2:C11"]
si_g = xf_cost_coeffs["SI"]["C2:C11"]
sk_g = xf_cost_coeffs["SK"]["C2:C11"]

albe_c = zeros(num_tech)
alde_c = zeros(num_tech)

alpha = vcat(albe_c, alde_c, at_a, be_a, cz_a, de_a, fr_a, hr_a, hu_a, nl_a, pl_a, ro_a, si_a, sk_a)

beta = vcat(albe_c, alde_c, at_b, be_b, cz_b, de_b, fr_b, hr_b, hu_b, nl_b, pl_b, ro_b, si_b, sk_b)

gamma = vcat(albe_c, alde_c, at_g, be_g, cz_g, de_g, fr_g, hr_g, hu_g, nl_g, pl_g, ro_g, si_g, sk_g)
