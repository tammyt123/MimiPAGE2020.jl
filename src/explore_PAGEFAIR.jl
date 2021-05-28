using Revise, MimiFAIR, MimiPAGE2020, Mimi, DataFrames

# include("../src/main_model_annual.jl")
# include("analysis/allscc/main_model_annual.jl") # run page on annual timestep

## load PAGE model
PAGE_m = MimiPAGE2020.get_model()
# m = getpage()
run(PAGE_m)
PAGE_mm = MimiPAGE2020.get_marginal_model(PAGE_m, year = 2030)
run(PAGE_mm)

## explore baseline outputs

# prtp and eta
PAGE_m[:EquityWeighting, :emuc_utilityconvexity] # eta = 1.66666....
PAGE_m[:EquityWeighting, :ptp_timepreference] # ptrp = 1.0333333....
PAGE_mm.base[:EquityWeighting, :emuc_utilityconvexity] # same as above
PAGE_mm.modified[:EquityWeighting, :emuc_utilityconvexity] # same as above
PAGE_mm.base[:EquityWeighting, :ptp_timepreference] # same as above
PAGE_mm.modified[:EquityWeighting, :ptp_timepreference] # same as above

# emissions
PAGE_mm[:CO2Cycle, :e_globalCO2emissions] # 0.1 in pulse year
PAGE_mm.base[:CO2Cycle, :e_globalCO2emissions]
PAGE_mm.modified[:CO2Cycle, :e_globalCO2emissions]
PAGE_mm.modified[:CO2Cycle, :e_globalCO2emissions] .- PAGE_mm.base[:CO2Cycle, :e_globalCO2emissions] # 7500 in pulse year

# temperature
PAGE_mm[:ClimateTemperature, :rt_g_globaltemperature] # temp delta only starts in the third time step? 0.0 for the second timestep (year = 2030)
PAGE_mm.base[:ClimateTemperature, :rt_g_globaltemperature]
PAGE_mm.modified[:ClimateTemperature, :rt_g_globaltemperature]
# PAGE_original_temp_delta = PAGE_mm.modified[:ClimateTemperature, :rt_g_globaltemperature] .- PAGE_mm.base[:ClimateTemperature, :rt_g_globaltemperature]

# total discounted impacts
PAGE_mm[:EquityWeighting, :td_totaldiscountedimpacts]
PAGE_mm.base[:EquityWeighting, :td_totaldiscountedimpacts]
PAGE_mm.modified[:EquityWeighting, :td_totaldiscountedimpacts]



PAGE_mm[:EquityWeighting, :td_totaldiscountedimpacts] / MimiPAGE2020.undiscount_scc(PAGE_mm.base, 2030)
PAGE_mm[:EquityWeighting, :td_totaldiscountedimpacts] / MimiPAGE2020.undiscount_scc(PAGE_mm.base, 2030)
MimiPAGE2020.compute_scc(PAGE_m, year = 2030)





# scc
MimiPAGE2020.compute_scc(PAGE_m, year = 2030, prtp = 0.03, eta = 0.) # 104.0209549634977
MimiPAGE2020.compute_scc(PAGE_m, year = 2030, prtp = 0.03, eta = 0., equity_weighting = false) # 104.0209549634977 ? 

MimiPAGE2020.compute_scc(year = 2030, prtp = 0.03, eta = 0.) # 104.0209549634977
MimiPAGE2020.compute_scc(year = 2030, prtp = 0.03, eta = 0., equity_weighting = false) # 104.0209549634977


## compare to PAGEFAIR

## load PAGE model
m = MimiPAGE2020.get_model()
# m = getpage()

eta = 0.0
prtp = 0.03

if eta !== nothing
    try
        set_param!(m, :emuc_utilityconvexity, eta)      # since eta is a default parameter in PAGE, we need to use `set_param!` if it hasn't been set yet
    catch e
        update_param!(m, :emuc_utilityconvexity, eta)   # or update_param! if it has been set
    end
end

if prtp !== nothing
    try
        set_param!(m, :ptp_timepreference, prtp * 100)      # since prtp is a default parameter in PAGE, we need to use `set_param!` if it hasn't been set yet
    catch e
        update_param!(m, :ptp_timepreference, prtp * 100)   # or update_param! if it has been set
    end
end

run(m)
# m[:EquityWeighting, :emuc_utilityconvexity] # 0.0
# m[:EquityWeighting, :ptp_timepreference] # 3.0

page_years = [2020, 2030, 2040, 2050, 2075, 2100, 2150, 2200, 2250, 2300]

## FAIR (modified version)
usg = "USG1"

FAIR = MimiFAIR.get_model(usg_scenario = usg)
run(FAIR)

FAIR[:temperature, :T] # annual temperature changes 1765-2300
fair_years = collect(1765:1:2300)

temperature = DataFrame(year = fair_years, T = FAIR[:temperature, :T])

input_temp = temperature[[year in page_years for year in fair_years], :T] # new baseline global average temperature
m[:ClimateTemperature, :rt_g_globaltemperature] # old baseline global average temperature
realized_temperature = input_temp * m[:ClimateTemperature, :ampf_amplification]' # new baseline regional realized temperatures
m[:ClimateTemperature, :rtl_realizedtemperature] # new baseline realized temperature

## set parameters in baseline PAGE model
set_param!(m, :rt_g_globaltemperature, input_temp)
set_param!(m, :rtl_realizedtemperature, realized_temperature)
run(m)

## create PAGE marginal model
pulse_size = 75000. # pulse size taken from PAGE2020
mm = Mimi.create_marginal_model(m, pulse_size) 
run(mm)

# note: PAGE adds a pulse to the pulse year of the size 75000 / period length, where period length = 0.5 * (year[i+1] - year[i-1]), where year[i] is the pulse year and years [i+1] and [i-1] correspond to the adjacent model years
# for year in page_years
#     println(75000/MimiPAGE2020.getperiodlength(year))
# end

## create FAIR marginal model
FAIR_mm = MimiFAIR.get_model(usg_scenario = usg)
run(FAIR_mm)

## set pulse year
pulse_year = 2030
pulse_year_index = findall((in)([pulse_year]), collect(1765:2300))
FAIR_pulse_size = pulse_size / MimiPAGE2020.getperiodlength(pulse_year) * 1/1000 * 12/44 # MtCO2 to GtC

## perturb CO2 emissions
new_emissions = FAIR_mm[:co2_cycle, :E_CO₂]
new_emissions[pulse_year_index] = new_emissions[pulse_year_index] .+ FAIR_pulse_size

MimiFAIR.update_param!(FAIR_mm, :E_CO₂, new_emissions)
run(FAIR_mm)
# new_temperature = FAIR_mm[:temperature, :T]

## input marginal FAIR temperature into marginal PAGE
new_temperature = DataFrame(year = fair_years, T = FAIR_mm[:temperature, :T])

new_input_temp = new_temperature[[year in page_years for year in fair_years], :T]
new_realized_temperature = new_input_temp * mm.base[:ClimateTemperature, :ampf_amplification]'

## set parameter in DICE model
MimiPAGE2020.set_param!(mm.modified, :rt_g_globaltemperature, new_input_temp)
MimiPAGE2020.set_param!(mm.modified, :rtl_realizedtemperature, new_realized_temperature)

# run(mm.modified)
run(mm)

scc = mm[:EquityWeighting, :td_totaldiscountedimpacts] / MimiPAGE2020.undiscount_scc(mm.base, pulse_year)
println(scc)

mm[:ClimateTemperature, :rt_g_globaltemperature] # 0.0, since it only updates intracomponent
mm.base[:ClimateTemperature, :rt_g_globaltemperature]
mm.modified[:ClimateTemperature, :rt_g_globaltemperature]

mm[:Discontinuity, :rt_g_globaltemperature] # around 10x smaller than original PAGE? check this
PAGE_mm[:Discontinuity, :rt_g_globaltemperature]

mm[:MarketDamages, :rtl_realizedtemperature]
PAGE_mm[:MarketDamages, :rtl_realizedtemperature]

mm[:EquityWeighting, :td_totaldiscountedimpacts]
PAGE_mm[:EquityWeighting, :td_totaldiscountedimpacts]

## plot

# global temperature delta
mm.base[:Discontinuity, :rt_g_globaltemperature]
mm.modified[:Discontinuity, :rt_g_globaltemperature]

mm.modified[:Discontinuity, :rt_g_globaltemperature] .- mm.base[:Discontinuity, :rt_g_globaltemperature]

plot(page_years, mm[:Discontinuity, :rt_g_globaltemperature]/1e9)
plot!(page_years, PAGE_mm[:Discontinuity, :rt_g_globaltemperature]/1e6)
# seems to be occuring one timestep later? 

# emissions
plot(fair_years, new_emissions*1000*44/12) # GtC to MtCO2
plot!(page_years, PAGE_mm.modified[:CO2Cycle, :e_globalCO2emissions])

# radiative forcing
plot(fair_years, (FAIR_mm[:co2_rf, :rf_co2] - FAIR[:co2_rf, :rf_co2])/1000*12/44, legend = :topleft)
plot!(page_years, PAGE_mm[:co2forcing, :f_CO2forcing])


