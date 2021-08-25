using Mimi
using CSV
using DataFrames

# specify model settings
function set_globalbools()
    global use_variability = false

    # set random seed to have similar variability development in the base and the marginal model.
    # set variability seed.
    if use_variability
        global varseed = rand(1:1000000000000)
    end

    global use_linear = false
    global use_logburke = false
    global use_logpopulation = false
    global use_logwherepossible = true
end

# set global values for technical configuration options
set_globalbools()

# get main_model file
include("analysis/allscc/main_model_annual.jl")
include("analysis/allscc/mcs_annual.jl")
include("analysis/allscc/compute_scc_annual.jl")



m = getpage()
m = get_model()
m = buildpage()

run(m)

m[:]