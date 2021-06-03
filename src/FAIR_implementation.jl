## function to get baseline PAGEFAIR model (i.e. input temperature vector from FAIR-NCEE into baseline PAGE2020)
function get_pagefair(;usg_scenario::String, prtp::Float64 = nothing, eta::Float64 = nothing)
    
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
    FAIR = MimiFAIR.get_model(usg_scenario = usg_scenario)
    run(FAIR)

    fair_years = collect(1765:1:2300)
    page_years = [2020, 2030, 2040, 2050, 2075, 2100, 2150, 2200, 2250, 2300] # this can be commented out when using the package
    temperature = DataFrame(year = fair_years, T = FAIR[:temperature, :T])

    input_temp = temperature[[year in page_years for year in fair_years], :T] # global average temperature input into baseline PAGE model
    realized_temperature = input_temp * m[:ClimateTemperature, :ampf_amplification]' # regional realized temperatures for baseline PAGE model: multiply input temp by regional amplification factors

    ## set parameters in baseline PAGE model
    set_param!(m, :rt_g_globaltemperature, input_temp)
    set_param!(m, :rtl_realizedtemperature, realized_temperature)
    run(m)

    return(m)

end

## function to get marginal DICE-FAIR model (i.e. input perturbed FAIR temperature vector into DICE2016)
## note: FAIR-NCEE must be loaded in environment first!!
function get_pagefair_marginal_model(;usg_scenario::String, pulse_year::Int, prtp::Float64 = nothing, eta::Float64 = nothing)
    
    ## create PAGE2020 marginal model
    m = MimiPAGE2020.get_pagefair(usg_scenario = usg_scenario, prtp = prtp, eta = eta)
    mm = Mimi.create_marginal_model(m, 1e9)
    run(mm)

    ## get perturbed FAIR temperature vector
    fair_years = collect(1765:1:2300)
    new_temperature = MimiFAIR.get_perturbed_fair_temperature(usg_scenario = usg_scenario, pulse_year = pulse_year)
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

## compute SCC from PAGEFAIR
function compute_scc_pagefair(;usg_scenario::String, pulse_year::Int, prtp::Float64, eta::Float64)
    mm = MimiPAGE2020.get_pagefair_marginal_model(usg_scenario = usg_scenario, pulse_year = pulse_year, prtp = prtp, eta = eta)
    scc = mm[:EquityWeighting, :td_totaldiscountedimpacts] / MimiPAGE2020.undiscount_scc(mm.base, pulse_year) * 1e6 # for 1 Gt FAIR pulse, since SCC is in millions
    return(scc)
end