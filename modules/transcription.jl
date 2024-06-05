using NLsolve
using NLopt
using LinearAlgebra
using ForwardDiff
using DifferentialEquations
using ComponentArrays
using Plots
using BenchmarkTools
using PolynomialRoots
using Distributions
using Measures
using InvertedIndices
using LaTeXStrings

abstract type TranscriptionModel end
include("./bulkmodels.jl")
include("./continuummodels.jl")
include("./modelselection.jl")
include("./rootsbackend.jl")


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

#Takes parameters and returns general case for analytic IC
function quasisteadyrate_analytic_general(k_i, K, α, β, θ, γ, Dtot, polymerasetot)
    p = polymerasetot ./Dtot
    cubicterm = α .*γ
    quadraticterm  = α .+ β ./Dtot .-γ .*(p .+θ ./ Dtot .+α)
    linearterm = -(p .+ α .+K ./ Dtot .-γ .* p)
    constantterm = p
    ic = polynomialroots(constantterm, linearterm, quadraticterm, cubicterm)
    rate = Dtot .*k_i .*ic .*(1 .-(β .*ic)/(K .-θ))/(1 .+γ .*ic)#Rate in nM RNA strands per hour
    return rate, ic 
end

#Takes same arguments as runtranscriptionmodel, returns quasisteady rate of RNA production
function quasisteadyrate_analytic(model::TranscriptionModel, params::AbstractArray{T1}, DNA, P) where {T1<:Number}
    (α, β, θ, γ) = (1, 0, 0, 0)
    K = (params.k_i+params.k_off)/params.k_on
    if typeof(model)<:Union{InitiationElongationModel, ContinuumModel}
        α = params.α
    end
    if typeof(model)<:ContinuumModel
        β = (α-1)*params.k_i*params.L/(params.k_on*params.N_all)
        θ = params.k_off/params.k_on
    end
    if typeof(model)<:LPTASEPModel
        γ = params.γ
    end
    (rate, ic) = quasisteadyrate_analytic_general(params.k_i, K, α, β, θ, γ, DNA, P)

    if typeof(model)<:BasicTASEPModel && ((β*ic)/(K-θ) .>0.5)[1]#Over approximation threshold
        rate = params.k_i*DNA*(K-θ)/(4*β)
    end
    if typeof(model)<:LPTASEPModel
        if 1+γ*(K-θ)/β > 0
            amax = (sqrt(1+γ*(K-θ)/β)-1)/(γ*(K-θ)/β)
        else
            amax = 0
        end
        if ((β*ic)/(K-θ).>amax)[1]#Over approximation threshold
        rate = params.k_i*DNA*((K-θ)/β)*amax*(1-amax)/(1+amax*γ*((K-θ)/β))
        end
    end
    return rate*params.N_all*1e-6#Output in mM NTP per hour
end

#Plot temporal outputs of model
function plottemporaloutput(model::TranscriptionModel, params, DNA, P, finaltime; kwargs...)
    plt = plot(xlabel = "Time (hrs)")
    plt = plottemporaloutput!(plt, model, params, DNA, P, finaltime; kwargs...)
    return plt
end

function plotQSrate_constantP(model::TranscriptionModel, params, DNArange, P; kwargs...)
    plt = plot(xlabel = "DNA (nM)", ylabel = "Transcription Rate (mM NTP/hr)")
    plt = plotQSrate_constantP!(plt, model, params, DNArange, P; kwargs...)
    return plt
end

function plotQSrate_constantD(model::TranscriptionModel, params, DNA, Prange; kwargs...)
    plt = plot(xlabel = "T7 RNA Polymerase per DNA", ylabel = "Transcription Rate per (mM NTP/ hr) per nM DNA")
    plt = plotQSrate_constantD!(plt, model, params, DNA, Prange; kwargs...)
    return plt
end

function plotQSrate_constantP!(plt, model::TranscriptionModel, params, DNArange, P; npoints = 20, modellabel = "", analytic=true)
    DNApoints = LinRange(DNArange[1],DNArange[2],npoints)
    QSrate = zeros(npoints)
    for (ind,DNA) in enumerate(DNApoints)
        if analytic
            QSrate[ind] = quasisteadyrate_analytic(model, params, DNA, P)
        else
            QSrate[ind] = dynamicquasisteadystate(model, params, DNA, P)
        end
    end
    plot!(plt,DNApoints,QSrate,label = modellabel)
end

function plotQSrate_constantD!(plt, model::TranscriptionModel, params, DNA, prange; npoints = 20, modellabel = "", analytic=true, color = :auto, normalizexaxis = true, normalizeyaxis = true)
    ppoints = LinRange(prange[1],prange[2],npoints)
    if normalizexaxis == true
        xscalingfactor = DNA
    else
        xscalingfactor = 1
    end

    if normalizeyaxis == true
        yscalingfactor = DNA
    else
        yscalingfactor = 1
    end

    QSrate = zeros(npoints)
    for (ind,p) in enumerate(ppoints)
        if analytic
            QSrate[ind] = quasisteadyrate_analytic(model, params, DNA, xscalingfactor*p)
        else
            QSrate[ind] = dynamicquasisteadystate(model, params, DNA, xscalingfactor*p)
        end
    end
    plot!(plt,ppoints,QSrate ./yscalingfactor,label = modellabel, linecolor = color)
end

#Functions for generating expressions common to differential models
function bulkprocessrates(param, DNA, P, IC, PE)
    r_on = param.k_on*DNA*P
    r_off = (param.k_off)*IC
    r_init = param.k_i*IC
    r_elong = (param.k_e/param.N_all)*PE
    return r_on, r_off, r_init, r_elong
end






