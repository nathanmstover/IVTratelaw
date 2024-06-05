abstract type PausingModel <: TranscriptionModel end
struct BasicPausingModel<:ContinuumModel end

NTPoutput(model::ContinuumModel, sol) = sol[end-1,:]

function plottemporaloutput!(plt, model::ContinuumModel, params, DNA, P, finaltime; modellabel = "")
    sol = runtranscriptionmodel(model, params, DNA, P, finaltime)
    NTPconsumed = NTPoutput(model,sol)
    RNAproduced = RNAoutput(model,sol)
    plot!(plt,sol.t,NTPconsumed, label = "NTP consumed: "*modellabel)
    plot!(plt,sol.t,RNAproduced, label = "RNA produced: "*modellabel)
    return plt
end

function continuumanimation(model::ContinuumModel, params, DNA, P, finaltime)
    sol = runtranscriptionmodel(model, params, DNA, P, finaltime)
    sol = sol(LinRange(0,finaltime,250))
    anim = @animate for i in 1:length(sol.t)
        plot(sol[i][4:end-2],ylims =(0,DNA))
    end
    gif(anim,"propagate.gif", fps = 50)
end

#Initialization of continuum (for basic TASEP model) initial conditions
function generateinitial(model::ContinuumModel, params::AbstractArray{T1},DNA,P) where {T1<:Number}
    nstatevars = Int(5+params.N_all)
    initial = zeros(T1,nstatevars)
    initial[1:nstatevars] = [DNA,P,zeros(nstatevars-2)...]
    return initial
end

#TASEP model
function modeldifferential!(model::BasicTASEPModel,du,u,param,t)
    N_all = Int(param.N_all)
    k_e_single = param.k_i*param.k_e

    DNA,P,IC = u[1:3]
    DNAtot = DNA+IC+sum(u[4:(end-1)])
    NTPconsumed = 0

    (r_on, r_off, r_initbulk) = bulkprocessrates(param, DNA, P, IC, 0)
    r_init = r_initbulk*(1-u[4]/(DNAtot))

    dDNA = -r_on+r_off+r_init
    dP = -r_on+r_off+k_e_single*u[N_all+3]
    dIC = r_on-r_off-r_init

    #DNA, P, and IC
    du[1:3] .= [dDNA, dP, dIC]
    #A1
    du[4] = r_init-k_e_single*u[4]*(1-u[5]/(DNAtot))
    #A2 through A_(N-1)
    for i in (5:N_all+2)
        added = k_e_single*u[i-1]*(1-u[i]/(DNAtot))
        lost = k_e_single*u[i]*(1-u[i+1]/(DNAtot))
        NTPconsumed += added
        du[i] = added-lost
    end
    #A_N
    added = k_e_single*u[N_all+2]*(1-u[N_all+3]/(DNAtot))
    lost = k_e_single*u[N_all+3]
    NTPconsumed += added
    du[N_all+3] = added-lost
    #NTP consumed
    du[N_all+4] =  NTPconsumed
    #Full RNA produced (NTP units)
    du[N_all+5] = N_all*lost
end

