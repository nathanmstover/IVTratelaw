function plotexperimentaldata!(plt,DPinputs,experimentaldata,experimentalnoise;
    colorscheme = :Dark2_5,
    extraDNA = [])

    nDNA = length(unique(DPinputs[:,1])) + length(extraDNA)
    uniquevals = reverse(unique(DPinputs[:,1]))
    for uniqueind in 1:length(uniquevals)
        color = palette(colorscheme,nDNA)[uniqueind]
        scatterx = DPinputs[DPinputs[:,1] .== uniquevals[uniqueind],2]
        scattery = experimentaldata[DPinputs[:,1] .== uniquevals[uniqueind]]
        noise = experimentalnoise[DPinputs[:,1] .== uniquevals[uniqueind]]
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


function analyzekineticdata(filename,N_all,N_plasmid,testmodelindex;meanK = 50, stdK = 20, massDNAunits = true, npoints = 50,showuncertainty = true, color = :Dark2_5, kwargs...)
    df = CSV.read(filename,DataFrame)
    (p_base,modellist,modelnamelist) = getbaseparams(N_all);
    Data = Matrix(df)
    extraDNA = []#For adding extra DNA lines to plot

    priorfunc = x->((x[2]-meanK)/stdK)^2
    priormatfull = [0 0 0 0 0;0 1/stdK^2 0 0 0;0 0 0 0 0;0 0 0 0 0;0 0 0 0 0]
    priormatlist = [priormatfull[1:2,1:2],priormatfull[1:3,1:3],priormatfull[1:4,1:4],priormatfull]
    priormat = priormatlist[testmodelindex]
    ###########You don't have to change anything below this line##################
    DPinputs = Data[:,1:2]
    experimentaldata = Data[:,3]
    experimentalnoise = Data[:,4]

    nDNA = length(unique(DPinputs[:,1])) + length(extraDNA)
    if massDNAunits
        DPinputs[:,1] .= DPinputs[:,1].* 1e9 ./(N_plasmid*607.4)#DNA concentration conversion
    else
        DPinputs[:,1] .= DPinputs[:,1]#DIFFERENT FROM MOST CONSTRUCTS AS UNITS ALREADY IN NM
    end
    DPinputs[:,2] .= DPinputs[:,2].*16#T7 concentration conversion

    testmodel = modellist[testmodelindex]
    testmodelname = modelnamelist[testmodelindex]

    normalizexaxis = false
    normalizeyaxis = false

    plt = plot()
    plotexperimentaldata!(plt,DPinputs,experimentaldata,experimentalnoise,extraDNA = extraDNA,colorscheme = color)
    (p_estimated,p_thintest,testcov) = multistartoptimizeparams(testmodel, experimentaldata, DPinputs, experimentalnoise, p_base; npoints = npoints, priorfunc = priorfunc, priormat = priormat, kwargs...)
    plotmodelprediction!(plt,testmodel,p_thintest,testcov,p_base,DPinputs,normalizexaxis = normalizexaxis,normalizeyaxis =normalizeyaxis, modelname = testmodelname, maximumT7 = 800, color = color, showuncertainty = showuncertainty, extraDNA = extraDNA)
    plot!(legend = :outerright,size = (450,550),leftmargin = 5mm,ylims = (0,150))
    display(plt)
    println("Parameter estimates (including 95% confidence intervals)")
    printuncertaintyestimates(p_thintest,testcov,N_all)
    pltelps2 = plotparametricellipse(p_thintest,testcov,title = testmodelname)
    display(pltelps2)
    return plt
end
