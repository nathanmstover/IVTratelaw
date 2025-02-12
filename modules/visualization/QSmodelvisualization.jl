function plotmodelprediction!(plt,model,params, covariancematrix, p_base, DPinputs;
    npoints = 35, 
    color = :blues,
    modelname="",
    normalizexaxis = true, 
    normalizeyaxis = true, 
    maximumT7 = 800,
    maximumrelativeT7 = 100, 
    extraDNA = [], 
    kwargs...)

    DNAs = reverse(vcat(unique(DPinputs[:,1]), extraDNA))
    colorpal = palette(color, length(DNAs))
    if normalizexaxis
        xlabel = "T7 RNA Polymerase per DNA"
        prange = (1e-5,maximumrelativeT7)
    else
        xlabel = "T7 RNA Polymerase (nM)"
        prange = (1e-5,maximumT7)
    end
    if normalizeyaxis
        ylabel = "Transcription Rate (mM NTP/ hr) per nM DNA"
    else
        ylabel = "Transcription Rate (mM NTP/ hr)"
    end
    plot!(plt, xlabel = xlabel, ylabel = ylabel, legendfontsize=9)
    for (ind,D) in enumerate(DNAs)
        plotQSrate_constantD!(plt, model, params, covariancematrix, p_base, D, prange; modellabel = "DNA = "*string(round(D, digits = 1))*" nM", npoints = npoints, color = colorpal[ind],normalizexaxis = normalizexaxis,normalizeyaxis = normalizeyaxis, kwargs...)
    end
    return plt
end

function plotQSrate_constantP(model::TranscriptionModel, params, DNArange, P; kwargs...)
    plt = plot(xlabel = "DNA (nM)", ylabel = "Transcription Rate (mM NTP/hr)")
    plt = plotQSrate_constantP!(plt, model, params, DNArange, P; kwargs...)
    return plt
end

function plotQSrate_constantD(model::TranscriptionModel, paramsls, covariancematrix, baseparams, DNA, Prange; kwargs...)
    plt = plot(xlabel = "T7 RNA Polymerase per DNA", ylabel = "Transcription Rate per (mM NTP/ hr) per nM DNA")
    plt = plotQSrate_constantD!(plt, model, paramsls, covariancematrix, baseparams, DNA, Prange; kwargs...)
    return plt
end

function plotQSrate_constantP!(plt, model::TranscriptionModel, params, DNArange, P; 
    npoints = 20, 
    modellabel = "", 
    analytic=true)

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

function plotQSrate_constantD!(plt, model::TranscriptionModel, paramsls, covariancematrix, baseparams, DNA, prange; 
    npoints = 20, 
    modellabel = "", 
    analytic=true, 
    color = :auto, 
    normalizexaxis = true, 
    normalizeyaxis = true, 
    showuncertainty = true, 
    thinparameters = true,
    alph = 0.05, 
    nmc = 10000)

    if thinparameters
        params = generatefullparameters(model,paramsls,baseparams)
    else
        params = paramsls
    end
    
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
    lower_pointwise_CB = zeros(npoints)
    upper_pointwise_CB = zeros(npoints)

    for (ind,p) in enumerate(ppoints)
        if analytic
            QSrate[ind] = quasisteadyrate_analytic(model, params, DNA, xscalingfactor*p)
        else
            QSrate[ind] = dynamicquasisteadystate(model, params, DNA, xscalingfactor*p)
        end
        if showuncertainty
            mcensemble = zeros(nmc)
            mean = paramsls
            d = MvNormal(mean, Hermitian(covariancematrix))
            for i in 1:nmc
                x = rand(d, 1)
                sampleparams = generatefullparameters(model,x,baseparams)
                mcensemble[i] = quasisteadyrate_analytic(model, sampleparams, DNA, xscalingfactor*p)
            end
            lower_pointwise_CB[ind] = percentile(mcensemble,100*alph/2)
            upper_pointwise_CB[ind] = percentile(mcensemble,100*(1-alph/2))
        else
            lower_pointwise_CB[ind] = 0
            upper_pointwise_CB[ind] = 0
        end
    end
    plot!(plt,ppoints,QSrate ./yscalingfactor,label = modellabel, linecolor = color,linewidth = 1.5)
    plot!(plt,ppoints, lower_pointwise_CB ./yscalingfactor, fillrange = upper_pointwise_CB ./yscalingfactor, fillalpha = 0.35,alpha=0.0, color = color,linewidth = 0.0,label="",z_order = :back)
end