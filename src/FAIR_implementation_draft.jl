using Revise, MimiFAIR, MimiPAGE2020, Mimi, DataFrames

## FIRST PASS OF FAIR IMPLEMENTATION (DO NOT USE FOR SCC -- see FAIR_compute_scc)

## load MimiPAGE2020
m = MimiPAGE2020.get_model()
run(m)

m[:ClimateTemperature, :rt_g_globaltemperature] # 10 element array
original_total_costs = m[:TotalCosts, :total_costs_peryear] # 10 x 8 array

page_years = [2020, 2030, 2040, 2050, 2075, 2100, 2150, 2200, 2250, 2300]

## FAIR (modified version)
FAIR = MimiFAIR.get_model(usg_scenario = "USG1")
run(FAIR)

FAIR[:temperature, :T] # global average temperature
fair_years = collect(1765:1:2300)

temperature = DataFrame(year = fair_years, T = FAIR[:temperature, :T])

input_temp = temperature[[year in page_years for year in fair_years], :T] # 10 element array

## calculate regional realized temperatures
m[:ClimateTemperature, :rtl_realizedtemperature]
m[:ClimateTemperature, :ampf_amplification]
m[:ClimateTemperature, :rt_g_globaltemperature] * m[:ClimateTemperature, :ampf_amplification]' == m[:ClimateTemperature, :rtl_realizedtemperature]
m[:ClimateTemperature, :rtl_realizedtemperature]

new_realized_temperature = input_temp * m[:ClimateTemperature, :ampf_amplification]'

## set new temperature parameters in PAGE model
# update_param!(m, :rt_g_globaltemperature, input_temp)
set_param!(m, :rt_g_globaltemperature, input_temp)
set_param!(m, :rtl_realizedtemperature, new_realized_temperature)
run(m)

## check and explore new outputs
m[:ClimateTemperature, :rt_g_globaltemperature] # original PAGE2020 temperature -- not updated
m[:Discontinuity, :rt_g_globaltemperature] # FAIR temperature -- updated
m[:SeaLevelRise, :rt_g_globaltemperature] # FAIR temperature -- updated

m[:ClimateTemperature, :rtl_realizedtemperature] # original PAGE temperature -- not updated
m[:NonMarketDamages, :rtl_realizedtemperature] # new realized temperature -- updated
# inter-component parameter connections have been updated, but inter-component parameters (i.e. within ClimateTemperature component) remain the original PAGE2020 values

new_total_costs = m[:TotalCosts, :total_costs_peryear]

new_total_costs - original_total_costs
sum(original_total_costs, dims = 2) # sum across regions for each time step
sum(new_total_costs, dims = 2) # sum across regions for each time step

sum(original_total_costs, dims = 2) .- sum(new_total_costs, dims = 2)
