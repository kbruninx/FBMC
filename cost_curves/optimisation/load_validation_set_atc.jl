using JuMP, HiGHS
#using Gurobi
using DataFrames, XLSX
using LinearAlgebra
using Statistics
using SparseArrays
using Formatting

function remove_missing(mx)
    mx[ismissing.(mx)] .= 0.0
    return mx
end

fbmc_zones =  ["ALBE", "ALDE", "AT", "BE", "CZ", "DE_LU", "FR", "HR", "HU", "NL", "PL", "RO", "SI", "SK"]
non_fbmc_zones =  ["CH", "GB", "ES", "IT_NORD"] # ["CH", "GB", "ES", "IT_NORD", "DK_1", "NO_2"]
all_zones = vcat(fbmc_zones, non_fbmc_zones)

borders = [
    ["FR", "CH"],
    ["FR", "ES"],
    ["FR", "GB"],
    ["FR", "IT_NORD"],
    ["AT", "IT_NORD"],
    ["AT", "CH"],
    ["BE", "GB"],
    ["NL", "GB"],
    ["SI", "IT_NORD"],
    #["NL", "DK_1"],
    #["NL", "NO_2"],
]

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
xf_generation_outages = XLSX.readxlsx("./data/validation/generation_outages.xlsx")
xf_ren_gen = XLSX.readxlsx("./data/validation/renewable_generation.xlsx")
xf_netpos = XLSX.readxlsx("./data/validation/net_positions.xlsx")
xf_ptdf = XLSX.readxlsx("./data/validation/ptdfs.xlsx")
xf_capacities = XLSX.readxlsx("./data/installed_capacities_corrected.xlsx")

xf_da_prices_non_fbmc = XLSX.readxlsx("./data/validation/day_ahead_prices_non_fbmc.xlsx")
xf_demand_non_fbmc = XLSX.readxlsx("./data/validation/demand_non_fbmc.xlsx")
xf_generation_non_fbmc = XLSX.readxlsx("./data/validation/generation_non_fbmc.xlsx")
xf_ren_gen_non_fbmc = XLSX.readxlsx("./data/validation/renewable_generation_non_fbmc.xlsx")
xf_capacities_non_fbmc = XLSX.readxlsx("./data/installed_capacities_non_fbmc_corrected.xlsx")
xf_generation_outages_non_fbmc = XLSX.readxlsx("./data/validation/generation_outages_non_fbmc.xlsx")


num_z = 14 # including ALBE and ALDE
num_z_non_fbmc = 4 #6
num_atc_border = 9 #11

num_val_t = 1346 # size of validation set
num_tech = 10
num_j = 133 # maximum amont on CNEs at a given time in the training dataset

coal_prices_g = vec(remove_missing(xf_fuel_prices["Sheet1"][sprintf1("B2:B%d", num_val_t+1)]))
oil_prices_g = vec(remove_missing(xf_fuel_prices["Sheet1"][sprintf1("C2:C%d", num_val_t+1)]))
gas_prices_g = vec(remove_missing(xf_fuel_prices["Sheet1"][sprintf1("D2:D%d", num_val_t+1)]))
eua_prices_g = vec(remove_missing(xf_fuel_prices["Sheet1"][sprintf1("E2:E%d", num_val_t+1)]))

coal_prices_g = convert(Vector{Float64}, coal_prices_g)
oil_prices_g = convert(Vector{Float64}, oil_prices_g)
gas_prices_g = convert(Vector{Float64}, gas_prices_g)
eua_prices_g = convert(Vector{Float64}, eua_prices_g)

at_gen_out_g = remove_missing(xf_generation_outages["AT"][sprintf1("B2:K%d", num_val_t+1)])
be_gen_out_g = remove_missing(xf_generation_outages["BE"][sprintf1("B2:K%d", num_val_t+1)])
cz_gen_out_g = remove_missing(xf_generation_outages["CZ"][sprintf1("B2:K%d", num_val_t+1)])
de_gen_out_g = remove_missing(xf_generation_outages["DE_LU"][sprintf1("B2:K%d", num_val_t+1)])
fr_gen_out_g = remove_missing(xf_generation_outages["FR"][sprintf1("B2:K%d", num_val_t+1)])
hr_gen_out_g = remove_missing(xf_generation_outages["HR"][sprintf1("B2:K%d", num_val_t+1)])
hu_gen_out_g = remove_missing(xf_generation_outages["HU"][sprintf1("B2:K%d", num_val_t+1)])
nl_gen_out_g = remove_missing(xf_generation_outages["NL"][sprintf1("B2:K%d", num_val_t+1)])
pl_gen_out_g = remove_missing(xf_generation_outages["PL"][sprintf1("B2:K%d", num_val_t+1)])
ro_gen_out_g = remove_missing(xf_generation_outages["RO"][sprintf1("B2:K%d", num_val_t+1)])
si_gen_out_g = remove_missing(xf_generation_outages["SI"][sprintf1("B2:K%d", num_val_t+1)])
sk_gen_out_g = remove_missing(xf_generation_outages["SK"][sprintf1("B2:K%d", num_val_t+1)])
ch_gen_out_g = remove_missing(xf_generation_outages_non_fbmc["CH"][sprintf1("B2:K%d", num_val_t+1)])
gb_gen_out_g = remove_missing(xf_generation_outages_non_fbmc["GB"][sprintf1("B2:K%d", num_val_t+1)])
es_gen_out_g = remove_missing(xf_generation_outages_non_fbmc["ES"][sprintf1("B2:K%d", num_val_t+1)])
it_gen_out_g = remove_missing(xf_generation_outages_non_fbmc["IT_NORD"][sprintf1("B2:K%d", num_val_t+1)])
dk_gen_out_g = remove_missing(xf_generation_outages_non_fbmc["DK_1"][sprintf1("B2:K%d", num_val_t+1)])
no_gen_out_g = remove_missing(xf_generation_outages_non_fbmc["NO_2"][sprintf1("B2:K%d", num_val_t+1)])

demand_g = remove_missing(xf_demand["Sheet1"][sprintf1("B2:O%d", num_val_t+1)]) # [z+t]
#demand_g_non_fbmc = remove_missing(xf_demand_non_fbmc["Sheet1"][sprintf1("B2:G%d", num_train_t+1)]) # [z+t]
demand_g_non_fbmc = remove_missing(xf_demand_non_fbmc["Sheet1"][sprintf1("B2:E%d", num_train_t+1)]) # [z+t]

ptdf_z_g = remove_missing(xf_ptdf["Sheet1"][sprintf1("D2:Q%d", num_val_t*num_j+1)])
ram_g = remove_missing(xf_ptdf["Sheet1"][sprintf1("C2:C%d", num_val_t*num_j+1)])

g_max_g = vec(xf_capacities["Sheet1"]["B2:O11"]) # [z+tech]
#g_max_g_non_fbmc = vec(xf_capacities_non_fbmc["Sheet1"]["B2:G11"]) # [z+tech]
g_max_g_non_fbmc = vec(xf_capacities_non_fbmc["Sheet1"]["B2:E11"]) # [z+tech]

ren_gen_g = remove_missing(xf_ren_gen["Sheet1"][sprintf1("B2:O%d", num_val_t+1)]) # [z+t]
#ren_gen_g_non_fbmc = remove_missing(xf_ren_gen_non_fbmc["Sheet1"][sprintf1("B2:G%d", num_train_t+1)]) # [z+t]
ren_gen_g_non_fbmc = remove_missing(xf_ren_gen_non_fbmc["Sheet1"][sprintf1("B2:E%d", num_train_t+1)]) # [z+t]

xf_cost_coeffs = XLSX.readxlsx("./cost_coefficients_atc_corrected.xlsx")

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
ch_a = xf_cost_coeffs["CH"]["A2:A11"]
gb_a = xf_cost_coeffs["GB"]["A2:A11"]
es_a = xf_cost_coeffs["ES"]["A2:A11"]
it_a = xf_cost_coeffs["IT_NORD"]["A2:A11"]
#dk_a = xf_cost_coeffs["DK_1"]["A2:A11"]
#no_a = xf_cost_coeffs["NO_2"]["A2:A11"]

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
ch_b = xf_cost_coeffs["CH"]["B2:B11"]
gb_b = xf_cost_coeffs["GB"]["B2:B11"]
es_b = xf_cost_coeffs["ES"]["B2:B11"]
it_b = xf_cost_coeffs["IT_NORD"]["B2:B11"]
#dk_b = xf_cost_coeffs["DK_1"]["B2:B11"]
#no_b = xf_cost_coeffs["NO_2"]["B2:B11"]

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
ch_g = xf_cost_coeffs["CH"]["C2:C11"]
gb_g = xf_cost_coeffs["GB"]["C2:C11"]
es_g = xf_cost_coeffs["ES"]["C2:C11"]
it_g = xf_cost_coeffs["IT_NORD"]["C2:C11"]
#dk_g = xf_cost_coeffs["DK_1"]["C2:C11"]
#no_g = xf_cost_coeffs["NO_2"]["C2:C11"]

albe_c = zeros(num_tech)
alde_c = zeros(num_tech)

alpha = vcat(albe_c, alde_c, at_a, be_a, cz_a, de_a, fr_a, hr_a, hu_a, nl_a, pl_a, ro_a, si_a, sk_a, ch_a, gb_a, es_a, it_a) #, dk_a, no_a)

beta = vcat(albe_c, alde_c, at_b, be_b, cz_b, de_b, fr_b, hr_b, hu_b, nl_b, pl_b, ro_b, si_b, sk_b, ch_b, gb_b, es_b, it_b) #, dk_b, no_b)

gamma = vcat(albe_c, alde_c, at_g, be_g, cz_g, de_g, fr_g, hr_g, hu_g, nl_g, pl_g, ro_g, si_g, sk_g, ch_g, gb_g, es_g, it_g) #, dk_g, no_g)
