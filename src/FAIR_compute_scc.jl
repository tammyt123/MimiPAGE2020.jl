using Revise, MimiFAIR, MimiPAGE2020, Mimi, DataFrames

## test: calculate SCC for original PAGE2020 without using compute_scc function

# # with original prtp and eta

# m = MimiPAGE2020.get_model()
# mm = MimiPAGE2020.get_marginal_model(year = 2030)
# scc_manual = mm[:EquityWeighting, :td_totaldiscountedimpacts] / MimiPAGE2020.undiscount_scc(mm.base, 2030)

# scc = MimiPAGE2020.compute_scc(m, year = 2030)

# # using constant 3% discount

# prtp = 0.03
# eta = 0.

# PAGE2020 = MimiPAGE2020.get_model()

# if eta !== nothing
#     try
#         set_param!(PAGE2020, :emuc_utilityconvexity, eta)      # since eta is a default parameter in PAGE, we need to use `set_param!` if it hasn't been set yet
#     catch e
#         update_param!(PAGE2020, :emuc_utilityconvexity, eta)   # or update_param! if it has been set
#     end
# end

# if prtp !== nothing
#     try
#         set_param!(PAGE2020, :ptp_timepreference, prtp * 100)      # since prtp is a default parameter in PAGE, we need to use `set_param!` if it hasn't been set yet
#     catch e
#         update_param!(PAGE2020, :ptp_timepreference, prtp * 100)   # or update_param! if it has been set
#     end
# end

# run(PAGE2020)
# PAGE2020_mm = MimiPAGE2020.get_marginal_model(PAGE2020, year = 2030)
# run(PAGE2020_mm)

# PAGE2020_mm[:EquityWeighting, :td_totaldiscountedimpacts] / MimiPAGE2020.undiscount_scc(PAGE2020_mm.base, 2030)
# MimiPAGE2020.compute_scc(PAGE2020, year = 2030, pulse_size = 75000., prtp = 0.03, eta = 0.)

# # results are approximately equal

####################################################
############### CALCULATE SCC ######################
####################################################

# 1. run baseline model (i.e. PAGE2020 + FAIR without emissions pulse)
# 2. run marginal model (i.e. PAGE2020 + FAIR with emissions pulse) (make sure to scale damages by pulse size)
# 3. calculate marginal damages and take PV to get SCC

# note: PAGE adds a pulse to the pulse year of the size 75000 / period length, where period length = 0.5 * (year[i+1] - year[i-1]), where year[i] is the pulse year and years [i+1] and [i-1] correspond to the adjacent model years
# i.e. pulse size equal to:
# for year in page_years
#     println(75000/MimiPAGE2020.getperiodlength(year))
# end

usg_scenario = ["USG1", "USG2", "USG3", "USG4", "USG5"]
usg = "USG1"

pulse_years = [2030, 2040, 2050, 2060]
pulse_year = 2030

eta = 0.0
prtp = 0.03

page_years = [2020, 2030, 2040, 2050, 2075, 2100, 2150, 2200, 2250, 2300]
fair_years = collect(1765:1:2300)

for usg in usg_scenario

    for pulse_year in pulse_years

        # --------------------------------------------------
        # run baseline model
        # --------------------------------------------------
        
        ## load PAGE model
        m = MimiPAGE2020.get_model()
        # m = getpage()

        # set prtp and eta
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

        ## FAIR (modified version)
        FAIR = MimiFAIR.get_model(usg_scenario = usg) # FAIR baseline model
        run(FAIR)

        # FAIR[:temperature, :T] # annual temperature changes 1765-2300        
        temperature = DataFrame(year = fair_years, T = FAIR[:temperature, :T])

        input_temp = temperature[[year in page_years for year in fair_years], :T] # global average temperature input into baseline PAGE model
        realized_temperature = input_temp * m[:ClimateTemperature, :ampf_amplification]' # regional realized temperatures for baseline PAGE model: multiply input temp by regional amplification factors

        ## set parameters in baseline PAGE model
        set_param!(m, :rt_g_globaltemperature, input_temp)
        set_param!(m, :rtl_realizedtemperature, realized_temperature)
        run(m)

        # --------------------------------------------------
        # run marginal model (i.e. add emissions pulse to FAIR)
        # --------------------------------------------------

        ## create PAGE marginal model
        # pulse_size = 75000. # pulse size in PAGE2020
        pulse_size = 1.0 # use 1 Mt pulse for now (PAGE2020 emissions are in Mt, so even though this is technically a "1 tonne pulse", it is computed as 1 Mt in PAGE)
        mm = Mimi.create_marginal_model(m, pulse_size) 
        run(mm)

        ## create FAIR marginal model
        FAIR_mm = MimiFAIR.get_model(usg_scenario = usg)
        run(FAIR_mm)

        ## set pulse year index
        pulse_year_index = findall((in)([pulse_year]), collect(1765:2300))

        ## set FAIR pulse size
        # FAIR_pulse_size = (pulse_size * 1/1000 * 12/44) / MimiPAGE2020.getperiodlength(pulse_year)  # MtCO2 to GtC
        # FAIR_pulse_size = (pulse_size * 1/1000 * 12/44)  # MtCO2 to GtC
        # FAIR_pulse_size = 1.0 * 12/44 # GtC to GtCO2 
        FAIR_pulse_size = 1/1000 * 12/44 # GtC to MtCO2 

        ## perturb FAIR CO2 emissions
        new_emissions = FAIR_mm[:co2_cycle, :E_CO₂]
        new_emissions[pulse_year_index] = new_emissions[pulse_year_index] .+ FAIR_pulse_size # add CO2 emissions to pulse year
        MimiFAIR.update_param!(FAIR_mm, :E_CO₂, new_emissions)
        run(FAIR_mm)
        # new_temperature = FAIR_mm[:temperature, :T]

        ## input perturbed FAIR temperature into marginal PAGE model
        new_temperature = DataFrame(year = fair_years, T = FAIR_mm[:temperature, :T])
        new_input_temp = new_temperature[[year in page_years for year in fair_years], :T]
        new_realized_temperature = new_input_temp * m[:ClimateTemperature, :ampf_amplification]' # amplification factors are the same for base and marginal models

        ## set parameters in marginal PAGE model
        MimiPAGE2020.set_param!(mm.modified, :rt_g_globaltemperature, new_input_temp)
        MimiPAGE2020.set_param!(mm.modified, :rtl_realizedtemperature, new_realized_temperature)

        # run(mm.modified)
        run(mm)

        # --------------------------------------------------
        # calculate SCC and print result
        # --------------------------------------------------

        # scc = mm[:EquityWeighting, :td_totaldiscountedimpacts] / MimiPAGE2020.undiscount_scc(mm.base, pulse_year) * 1e6 / 1e9 # millions to dollars, Gt to t (for 1 Gt FAIR pulse)
        scc = mm[:EquityWeighting, :td_totaldiscountedimpacts] / MimiPAGE2020.undiscount_scc(mm.base, pulse_year) # for 1 Mt FAIR pulse
        println(scc)

    end
end
