abstract type TranscriptionModel end

#Function for generating expressions common to differential models
function bulkprocessrates(param, DNA, P, IC, PE)
    r_on = param.k_on*DNA*P
    r_off = (param.k_off)*IC
    r_init = param.k_i*IC
    r_elong = (param.k_e/param.N_all)*PE
    return r_on, r_off, r_init, r_elong
end

include("./bulkmodels.jl")
include("./continuummodels.jl")
include("./analyticapproximations.jl")
include("./dynamicsimulation.jl")





