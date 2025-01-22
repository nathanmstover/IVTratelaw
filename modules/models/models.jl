abstract type TranscriptionModel end

#Function for generating expressions common to differential models
function bulkprocessrates(param, DNA, P, IC, PE)
    r_on = param.k_on*DNA*P
    r_off = (param.k_off)*IC
    r_init = param.k_i*IC
    r_elong = (param.k_e/param.N_all)*PE
    return r_on, r_off, r_init, r_elong
end

function getbaseparams(N_all)
    #Basic parameters (from Stover et al, 2023)
    k_on = 10^2.30 #P-DNA binding (hours^-1 nM^-1)
    k_off = 10^3.74 #P-DNA binding (inverse hours)
    k_i = 10^3.61 #Initiation rate (inverse hours)
    k_e = 10^5.20 #Elongation rate (inverse hours)

    #Basic pausing parameters (Klumpp and Hwa, 2008)
    τ = 1/3600 #Pausing time (hours)
    f = 0.1*3600 #Pausing frequency (inverse hours)

    #Basic T7 RNA polymerase size parameter (My best guess + Klumpp and Hwa, 2008)
    L = 25
    p_base = ComponentVector((k_on = k_on, k_off = k_off, k_i = k_i, k_e = k_e, L = L, f = f, τ = τ, N_all = N_all))

    initiationlimited = InitiationLimitedModel()
    initiationelongation = InitiationElongationModel()
    basicTASEP = BasicTASEPModel()
    longTASEP = LPTASEPModel();
    modellist = [initiationlimited,initiationelongation,basicTASEP,longTASEP]
    modelnamelist = ["Initiation limited","Initiation-elongation","Basic TASEP","Long Pause TASEP"]

    return p_base, modellist, modelnamelist
end

include("./bulkmodels.jl")
include("./continuummodels.jl")
include("./analyticapproximations.jl")
include("./dynamicsimulation.jl")




