abstract type BulkModel <: TranscriptionModel end
struct InitiationLimitedModel<:BulkModel end
struct InitiationElongationModel<:BulkModel end

function plottemporaloutput!(plt, model::BulkModel, params, DNA, P, finaltime; modellabel = "")
    sol = dynamictranscriptionmodel(model, params, DNA, P, finaltime)
    NTPconsumed = NTPoutput(model,sol)
    plot!(plt,sol.t,NTPconsumed, label = "NTP consumed: "*modellabel)
    return plt
end

#Functions for generating initial conditions
function generateinitial(model::InitiationLimitedModel, params::AbstractArray{T1},DNA,P) where {T1<:Number}
    nstatevars = 4
    initial = zeros(T1,nstatevars)
    initial[1:nstatevars] = [DNA,P,0,0]
    return initial
end

function generateinitial(model::InitiationElongationModel, params::AbstractArray{T1},DNA,P) where {T1<:Number}
    nstatevars = 5
    initial = zeros(T1,nstatevars)
    initial[1:nstatevars] = [DNA,P,0,0,0]
    return initial
end

function modeldifferential!(model::InitiationLimitedModel,du,u,param,t)
    DNA,P,IC,RNA = u
    (r_on, r_off, r_init) = bulkprocessrates(param, DNA, P, IC, 0)
    dDNA = dP = -r_on+r_off+r_init
    dIC = r_on-r_off-r_init
    #DNA, Polymerase, and IC concentrations
    du[1:3] .= [dDNA,dP,dIC]
    #NTP consumed and RNA produced
    du[4] = param.N_all*r_init
end

function modeldifferential!(model::InitiationElongationModel,du,u,param,t)
    DNA,P,IC,PE,RNA = u
    (r_on, r_off, r_init, r_elong) = bulkprocessrates(param, DNA, P, IC, PE)
    dDNA = -r_on+r_off+r_init
    dP = -r_on+r_off+r_elong
    dIC = r_on-r_off-r_init
    dPE = r_init-r_elong
    #DNA, Polymerase, IC, and PE concentrations
    du[1:4] .= [dDNA, dP, dIC, dPE]
    #NTP consumed and RNA produced
    du[5] =  param.N_all*r_elong
end

