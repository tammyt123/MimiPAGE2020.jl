using Mimi

include("analysis/allscc/main_model_annual.jl")

m = getpage()
run(m)s
explore(m)


## to run PAGE2020 on an annual timestep
# using Mimi
include("../src/main_model_annual.jl")
# include("../src/main_model_annual.jl")

# m = getpage()
# run(m)
# explore(m)