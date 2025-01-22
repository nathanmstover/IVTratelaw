function modelprediction(model,estimatedparams, DPpoints, baseparams)
    fullparams = generatefullparameters(model,estimatedparams,baseparams)
    generatesyntheticdata(model,fullparams,DPpoints,zeros(size(DPpoints)[1]))
end

function generatesyntheticdata(model,params,DPinputs,noise)
    ndatapoints = size(DPinputs)[1]
    syntheticdata = []
    for ind in 1:ndatapoints
        (DNA,T7) = DPinputs[ind,:]
        syntheticpoint = quasisteadyrate_analytic(model, params, DNA, T7) + noise[ind]*randn()
        append!(syntheticdata,syntheticpoint)
    end
    return syntheticdata
end

function generatefullparameters(model::InitiationLimitedModel,estimatedparams,baseparams)
    return ComponentVector((k_i = estimatedparams[1], K = estimatedparams[2], N_all = baseparams.N_all)) 
end
function generatefullparameters(model::T,estimatedparams,baseparams) where T<:InitiationElongationModel
    k_i = estimatedparams[1]
    β = baseparams.L*k_i^2/(baseparams.k_e*baseparams.k_on)
    θ = baseparams.k_off/baseparams.k_on
    return ComponentVector((k_i = k_i, K = estimatedparams[2], α = estimatedparams[3], β = β, θ = θ, N_all = baseparams.N_all)) 
end
function generatefullparameters(model::T,estimatedparams,baseparams) where T<:BasicTASEPModel
    k_i = estimatedparams[1]
    β = baseparams.L*k_i^2/(baseparams.k_e*baseparams.k_on)
    θ = baseparams.k_off/baseparams.k_on
    return ComponentVector((k_i = k_i, K = estimatedparams[2], α = estimatedparams[3], β = estimatedparams[4], θ = θ, N_all = baseparams.N_all)) 
end
function generatefullparameters(model::LPTASEPModel,estimatedparams,baseparams)
    k_i = estimatedparams[1]
    β = baseparams.L*k_i^2/(baseparams.k_e*baseparams.k_on)
    θ = baseparams.k_off/baseparams.k_on
    return ComponentVector((k_i = k_i, K = estimatedparams[2], α = estimatedparams[3], β = estimatedparams[4], γ = estimatedparams[5], θ = θ, N_all = baseparams.N_all)) 
end

#Take list of DNA inputs and list of polymerase inputs, outputs a list of DNA,polymerase input pairs for a full grid experimental design
function generategridinputs(DNApoints,Ppoints)
    (hcat([hcat([[D,P] for P in Ppoints]...) for D in DNApoints]...)')
end



