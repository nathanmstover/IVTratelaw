abstract type ContinuumModel <: TranscriptionModel end
struct BasicTASEPModel<:ContinuumModel end
struct LPTASEPModel<:ContinuumModel end

NTPoutput(model::ContinuumModel, sol) = sol[end-1,:] .* 1e-6 #RNA in mM NTP equivalents

function plottemporaloutput!(plt, model::ContinuumModel, params, DNA, P, finaltime; modellabel = "")
    sol = dynamictranscriptionmodel(model, params, DNA, P, finaltime)
    NTPconsumed = NTPoutput(model,sol)
    RNAproduced = RNAoutput(model,sol)
    plot!(plt,sol.t,NTPconsumed, label = "NTP consumed: "*modellabel)
    plot!(plt,sol.t,RNAproduced, label = "RNA produced: "*modellabel)
    return plt
end

function continuumanimation(model::ContinuumModel, params, DNA, P, finaltime)
    sol = dynamictranscriptionmodel(model, params, DNA, P, finaltime)
    sol = sol(LinRange(0,finaltime,250))
    anim = @animate for i in 1:length(sol.t)
        plot(sol[i][4:end-2] ./DNA,ylims =(0,1.1),xlabel = "Sequence location index", ylabel = "Fraction of sites occupied", label = "")
    end
    gif(anim,"figures/propagate.gif", fps = 50)
end

function densitydistributionfinal(model::ContinuumModel, params, DNA, P; kwargs...)
    plt = plot(xlabel = "DNA position", ylabel = "Polymerase occupancy fraction")
    densitydistributionfinal!(plt, model, params, DNA, P; kwargs...)
    return plt
end

function densitydistributionfinal!(plt, model::ContinuumModel, params, DNA, P; plotlabel = "")
    sol = dynamictranscriptionmodel(model, params, DNA, P, 1)#Assuming 1 hour is enough for SS
    plot!(sol[end][4:end-2] ./DNA, ylims = (0,1.1), label = plotlabel,linewidth = 2)
    return plt
end

#Initialization of continuum (for basic TASEP model) initial conditions
function generateinitial(model::ContinuumModel, params::AbstractArray{T1},DNA,P) where {T1<:Number}
    M = Int(round(params.N_all/params.L))
    nstatevars = Int(5+M)
    initial = zeros(T1,nstatevars)
    initial[1:nstatevars] = [DNA,P,zeros(nstatevars-2)...]
    return initial
end

function blockingfactor(model::BasicTASEPModel,segmentconcentration,DNAtot,param)
    a = segmentconcentration/(DNAtot)
    return (1-a)
end

function blockingfactor(model::LPTASEPModel,segmentconcentration,DNAtot,param)
    k_e = param.k_e/param.L
    a = segmentconcentration/(DNAtot)
    return (1-a)/(1+a*k_e*param.f*param.τ^2/(1+f*τ))
end

#Basic TASEP model
function modeldifferential!(model::ContinuumModel,du,u,param,t)
    M = Int(round(param.N_all/param.L))
    k_e = param.k_e/param.L
    DNA,P,IC = u[1:3]
    DNAtot = DNA+IC
    NTPconsumed = 0

    (r_on, r_off, r_initbulk) = bulkprocessrates(param, DNA, P, IC, 0)

    r_init = r_initbulk*blockingfactor(model,u[4],DNAtot,param)

    dDNA = -r_on+r_off+r_init
    dP = -r_on+r_off+k_e*u[M+3]
    dIC = r_on-r_off-r_init

    #DNA, P, and IC
    du[1:3] .= [dDNA, dP, dIC]
    #A1
    du[4] = r_init-k_e*u[4]*blockingfactor(model,u[5],DNAtot,param)
    #A2 through A_(N-1)
    for i in (5:M+2)
        added = k_e*u[i-1]*blockingfactor(model,u[i],DNAtot,param)
        lost = k_e*u[i]*blockingfactor(model,u[i+1],DNAtot,param)
        NTPconsumed += param.L*added
        du[i] = added-lost
    end
    #A_N
    added = k_e*u[M+2]*blockingfactor(model,u[M+3],DNAtot,param)
    lost = k_e*u[M+3]
    NTPconsumed += param.L*added
    du[M+3] = added-lost
    #NTP consumed
    du[M+4] =  NTPconsumed
    #Full RNA produced (NTP units)
    du[M+5] = param.N_all*lost
end