#Defining Functions for Observables
RNAoutput(model::TranscriptionModel, sol) = sol[end,:] .* 1e-6 #RNA in mM NTP equivalents
NTPoutput(model::TranscriptionModel, sol) = sol[end,:] .* 1e-6 #NTP in mM

function dynamictranscriptionmodel(model::TranscriptionModel, params::AbstractArray{T1}, DNA, P, finaltime) where {T1<:Number}
    initial = generateinitial(model,params,DNA,P)
    differentialeqn = (du,u,param,t) -> modeldifferential!(model,du,u,param,t)
    f = ODEFunction(differentialeqn)
    prob = ODEProblem(f,initial,(0,finaltime),params)
    sol = solve(prob,Tsit5())
    return sol
end

function dynamicquasisteadystate(model::TranscriptionModel, params::AbstractArray{T1}, DNA, P) where {T1<:Number}
    #Assuming that 1 hour is enough to reach QSS
    10*(NTPoutput(model,dynamictranscriptionmodel(model, params, DNA, P, 0.2))[end] - NTPoutput(model,dynamictranscriptionmodel(model, params, DNA, P, 0.1))[end])
end