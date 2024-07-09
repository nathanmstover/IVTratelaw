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
