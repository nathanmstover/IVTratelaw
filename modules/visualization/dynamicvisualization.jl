#Plot temporal outputs of model
function plottemporaloutput(model::TranscriptionModel, params, DNA, P, finaltime; kwargs...)
    plt = plot(xlabel = "Time (hrs)")
    plt = plottemporaloutput!(plt, model, params, DNA, P, finaltime; kwargs...)
    return plt
end

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

"""
    relative(f, r; sp)

Generate relative coordinates for plotting.""" 
function relative(f, r; sp)
    p = plot!()
    lims = f(p[sp])
    return lims[1] + r * (lims[2]-lims[1])
end
relativex(r; sp::Int=1) = relative(Plots.xlims, r; sp=sp)
relativey(r; sp::Int=1) = relative(Plots.ylims, r; sp=sp)