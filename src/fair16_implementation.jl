#######################################################################################################################
# LOAD BASELINE PAGEFAIR MODEL
########################################################################################################################
# Description: Return PAGE2020 model with temperature vector set to baseline FAIR model. Note: in order to run this, FAIR-NCEE must
#              be loaded in the environment.
#
# Function Arguments:
#
#       ar6_scenario:     ar6 scenario
#       prtp:             Pure rate of time preference, for calculating discount rate.
#       eta:              Eta parameter, for calculating discount rate.
#----------------------------------------------------------------------------------------------------------------------

function get_pagefair(;ar6_scenario::String="ssp245", prtp = nothing, eta = nothing)
    
    ## load baseline PAGE model
    m = MimiPAGE2020.get_model()
    run(m)

    # set prtp and eta in base model
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

    ## load baseline FAIR-NCEE
    FAIR = MimiFAIRv1_6_2.get_model(ar6_scenario = ar6_scenario)
    run(FAIR)

    fair_years = collect(1750:1:2300)
    page_years = [2020, 2030, 2040, 2050, 2075, 2100, 2150, 2200, 2250, 2300] # this can be commented out when using the package
    temperature = DataFrame(year = fair_years, T = FAIR[:temperature, :T])

    input_temp = temperature[[year in page_years for year in fair_years], :T] # global average temperature input into baseline PAGE model
    realized_temperature = input_temp * m[:ClimateTemperature, :ampf_amplification]' # regional realized temperatures for baseline PAGE model: multiply input temp by regional amplification factors

    ## set parameters in baseline PAGE model
    MimiPAGE2020.set_param!(m, :rt_g_globaltemperature, input_temp)
    MimiPAGE2020.set_param!(m, :rtl_realizedtemperature, realized_temperature)
    run(m)

    return(m)

end

#######################################################################################################################
# GET MARGINAL PAGEFAIR MODEL
########################################################################################################################
# Description: Return marginal (perturbed) PAGE2020 model with temperature vector set to perturbed FAIR model. 
#              Note: in order to run this, FAIR-NCEE must be loaded in the environment.
#
# Function Arguments:
#
#       ar6_scenario:     ar6 scenario
#       prtp:             Pure rate of time preference, for calculating discount rate.
#       eta:              Eta parameter, for calculating discount rate.
#       gas:              Gas to perturb (:CO2, :CH4, or :N2O).
#       pulse_year:       Pulse year (for SC-GHG calculation).
#       pulse_size:       Pulse size (defaults to 1.0).
#----------------------------------------------------------------------------------------------------------------------

function get_pagefair_marginal_model(;ar6_scenario::String="ssp245", pulse_year::Int, prtp = nothing, eta = nothing, gas::Symbol=:CO2, pulse_size::Float64=1.0)
    
    ## create PAGE2020 marginal model
    m = MimiPAGE2020.get_pagefair(ar6_scenario = ar6_scenario, prtp = prtp, eta = eta)
    # mm = Mimi.create_marginal_model(m, (pulse_size * 1e9)) # multiply by 1e9 since FAIR units are Gt
    mm = Mimi.create_marginal_model(m, pulse_size)
    run(mm)

    ## get perturbed FAIR temperature vector
    fair_years = collect(1750:1:2300)
    new_temperature = MimiFAIRv1_6_2.get_perturbed_fair_temperature(ar6_scenario = ar6_scenario, pulse_year = pulse_year, pulse_size = pulse_size, gas = gas)
    new_temperature_df = DataFrame(year = fair_years, T = new_temperature)

    ## input perturbed FAIR temperature into marginal PAGE model
    new_input_temp = new_temperature_df[[year in page_years for year in fair_years], :T]
    new_realized_temperature = new_input_temp * m[:ClimateTemperature, :ampf_amplification]' # amplification factors are the same for base and marginal models

    ## set parameters in marginal PAGE model
    MimiPAGE2020.set_param!(mm.modified, :rt_g_globaltemperature, new_input_temp)
    MimiPAGE2020.set_param!(mm.modified, :rtl_realizedtemperature, new_realized_temperature)

    run(mm)

    return(mm)

end


#######################################################################################################################
# COMPUTE SC-GHG
########################################################################################################################
# Description: Compute SC-GHGs. Note: in order to run this, FAIR-NCEE must be loaded in the environment.
#
# Function Arguments:
#
#       ar6_scenario:     ar6 scenario
#       prtp:             Pure rate of time preference, for calculating discount rate.
#       eta:              Eta parameter, for calculating discount rate.
#       gas:              Gas to perturb (:CO2, :CH4, or :N2O).
#       pulse_year:       Pulse year (for SC-GHG calculation).
#       pulse_size:       Pulse size (defaults to 1.0).
#----------------------------------------------------------------------------------------------------------------------

function compute_scghg_pagefair(;ar6_scenario::String="ssp245", pulse_year::Int, prtp::Float64, eta::Float64, pulse_size::Float64=1.0, gas::Symbol=:CO2)

    mm = MimiPAGE2020.get_pagefair_marginal_model(ar6_scenario = ar6_scenario, pulse_year = pulse_year, prtp = prtp, eta = eta, gas = gas, pulse_size = pulse_size)

    if gas == :CO2
        scghg = mm[:EquityWeighting, :td_totaldiscountedimpacts] / MimiPAGE2020.undiscount_scc(mm.base, pulse_year) / 1e3 * 12/44 # fair1.6 pulse is 1GtC, scc is in millions. need to convert to CO2 and divide by 1e3
    else
        scghg = mm[:EquityWeighting, :td_totaldiscountedimpacts] / MimiPAGE2020.undiscount_scc(mm.base, pulse_year)
    end

    return(scghg)

end
