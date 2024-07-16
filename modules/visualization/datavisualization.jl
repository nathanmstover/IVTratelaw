function plotexperimentaldata!(plt,DPinputs,experimentaldata,experimentalnoise;
    colorscheme = :Dark2_5,
    extraDNA = [])

    nDNA = length(unique(DPinputs[:,1])) + length(extraDNA)

    for uniqueind in 1:length(unique(DPinputs[:,1]))
        color = palette(colorscheme,nDNA)[uniqueind]
        scatterx = DPinputs[DPinputs[:,1] .== unique(DPinputs[:,1])[uniqueind],2]
        scattery = experimentaldata[DPinputs[:,1] .== unique(DPinputs[:,1])[uniqueind]]
        noise = experimentalnoise[DPinputs[:,1] .== unique(DPinputs[:,1])[uniqueind]]
        plt = scatter!(scatterx,scattery, label = "",color = color, yerror = noise)
    end
end

function generateplotsyntheticdata(model,params,DPinputs; noise = zeros(size(DPinputs[1])), kwargs...)
    plt = plot()
    #plotmodelprediction!(plt,model,params,DPinputs,zeros(2,2), p_base; color = :greens, showuncertainty = false, kwargs...)
    syntheticdata = generatesyntheticdata(model,params,DPinputs,noise)
    plotsyndata_constantD!(plt,syntheticdata, DPinputs, noise; kwargs...)
    return (plt,syntheticdata,DPinputs)
end

function plotsyndata_constantD!(plt,syntheticdata, DPinputs, noise; 
    color = :greens, 
    modelname="", 
    normalizexaxis = true, 
    normalizeyaxis = true,
    extraDNA = [])

    DNApoints = unique(DPinputs[:,1])
    colorpal = palette(color, length(DNApoints)+length(extraDNA))
    for (DNAind,DNA) in enumerate(DNApoints)
        T7points = DPinputs[isapprox.(DPinputs[:,1],DNA),2]
        transcriptionrate = syntheticdata[isapprox.(DPinputs[:,1],DNA)]
        if normalizexaxis
            xaxis = T7points ./ DNA
        else
            xaxis = T7points
        end
        if normalizeyaxis
            yaxis = transcriptionrate ./ DNA
            plotnoise = noise/DNA
        else
            yaxis = transcriptionrate
            plotnoise = noise
        end
        scatter!(plt,xaxis, yaxis, yerror = plotnoise, mc = colorpal[DNAind], label = "")
    end
end

function plotonboxconstraints(DPinputs;  
    upperDNA = 28.5, #nM
    lowerDNA = 1.0, #nM
    upperP = 384.0, #nM
    lowerP = 40.0, #nM
    DNAconcentration = 504,#stock solution nM
    T7concentration = 3200,#stock solution nM
    maxvolumeratio = 0.28)#fraction of vial
        
    maxvolumefunction = (DNAadded) -> T7concentration*(maxvolumeratio - DNAadded/DNAconcentration)
    DNArange = LinRange(lowerDNA,upperDNA,100)
    scatter(DPinputs[:,1],DPinputs[:,2], label = "", xlabel = "Added DNA (nM)", ylabel = "Added RNA polymerase (nM)", legend = :outertop, size = (500,500))
    vline!([lowerDNA,upperDNA], label = "DNA box constraints")
    hline!([lowerP,upperP], label = "RNA Polymerase box constraints")
    plot!(DNArange,[maxvolumefunction(x) for x in DNArange], label = "Maximum volume constraint")
    plot!(DNArange,lowerP*ones(length(DNArange)), fillrange = [max(lowerP,min(maxvolumefunction(x),upperP)) for x in DNArange], fillalpha = 0.35,alpha=0.0,linewidth = 0.0,label="Feasible Region",z_order = :back)

end
