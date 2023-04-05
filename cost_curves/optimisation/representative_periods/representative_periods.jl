using RepresentativePeriodsFinder
using JuMP
using Gurobi

config_file = "./representative_periods/representative_periods_config.yaml"

pf = PeriodsFinder(config_file)

optimizer = optimizer_with_attributes(Gurobi.Optimizer, "TimeLimit" => 300) # seconds specifies time out limit
pf = find_representative_periods(pf, optimizer=optimizer)