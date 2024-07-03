function generatesyntheticdata(model,params,DPinputs,noise)
    ndatapoints = size(DPinputs)[1]
    syntheticdata = []
    for ind in 1:ndatapoints
        (DNA,T7) = DPinputs[ind,:]
        syntheticpoint = quasisteadyrate_analytic(model, params, DNA, T7) + noise*randn()
        append!(syntheticdata,syntheticpoint)
    end
    return syntheticdata
end

function modelprediction(model,estimatedparams, DPpoints, baseparams)
    fullparams = generatefullparameters(model,estimatedparams,baseparams)
    generatesyntheticdata(model,fullparams,DPpoints,0)
end

function getsensitivitymatrix(model,estimatedparams, DPpoints, baseparams)
    fun = p->modelprediction(model,p, DPpoints, baseparams)
    ForwardDiff.jacobian(fun,estimatedparams)'
end

function getexperimentalvariancematrix(model,estimatedparams, DPpoints, baseparams, noise; homoscedastic = true)
    ndatapoints = size(DPpoints)[1]
    if homoscedastic
        return noise^2*I(ndatapoints)
    end
end

function getDoptimalitymeasure(model,estimatedparams, DPpointsflattened, baseparams, noise; priormatrix = zeros(length(estimatedparams),length(estimatedparams)), homoscedastic = true)
    nexperiments = Int(length(DPpointsflattened)/2)
    DPpoints = reshape(DPpointsflattened,(nexperiments,2))
    M = getinformationmatrix(model,estimatedparams, DPpoints, baseparams, noise, priormatrix = priormatrix, homoscedastic = homoscedastic)
    return det(inv(M))
end



function getinformationmatrix(model,estimatedparams, DPpoints, baseparams, noise; priormatrix = zeros(length(estimatedparams),length(estimatedparams)), homoscedastic = true)
    V = getexperimentalvariancematrix(model,estimatedparams, DPpoints, baseparams, noise; homoscedastic = homoscedastic)
    S = getsensitivitymatrix(model,estimatedparams, DPpoints, baseparams)
    return S*inv(V)*S' .+ priormatrix
end

function getcovariancematrix(model,estimatedparams, DPpoints, baseparams, noise; priormatrix = homoscedastic = true)
    M = getinformationmatrix(model,estimatedparams, DPpoints, baseparams, noise, priormatrix = 0, homoscedastic = true)
    return inv(M)
end

function generatefullparameters(model::InitiationLimitedModel,estimatedparams,baseparams)
    return ComponentVector((k_i = estimatedparams[1], L = baseparams.L, N_all = baseparams.N_all, k_off = baseparams.k_off, k_on = baseparams.k_on)) 
end
function generatefullparameters(model::T,estimatedparams,baseparams) where T<:Union{InitiationElongationModel,BasicTASEPModel}
    return ComponentVector((k_i = estimatedparams[1], α = estimatedparams[2], L = baseparams.L, N_all = baseparams.N_all, k_off = baseparams.k_off, k_on = baseparams.k_on)) 
end
function generatefullparameters(model::LPTASEPModel,estimatedparams,baseparams)
    return ComponentVector((k_i = estimatedparams[1], α = estimatedparams[2], γ = estimatedparams[3], L = baseparams.L, N_all = baseparams.N_all, k_off = baseparams.k_off, k_on = baseparams.k_on)) 
end

function modelresidual(model,estimatedparams, syntheticdata, DPpoints, noise, baseparams)
    predicteddata = modelprediction(model,estimatedparams, DPpoints, baseparams)
    norm(predicteddata .- syntheticdata)/noise
end

function fittingobjective(x::Vector, grad::Vector, model, syntheticdata, DPpoints, noise, baseparams)
    if length(grad) > 0
        fun = p->modelresidual(model, p, syntheticdata, DPpoints, noise, baseparams)
        gnew = ForwardDiff.gradient(fun,x)
        for i in 1:length(grad)
            grad[i] = gnew[i]
        end
    end
    residual = modelresidual(model, x, syntheticdata, DPpoints, noise, baseparams)
    return residual
end

function localoptimizeparams(model, syntheticdata, DPpoints, noise, baseparams, initguess, nparams, upperbounds, lowerbounds)
    opt = Opt(:LD_SLSQP, nparams)
    opt.lower_bounds = lowerbounds
    opt.upper_bounds = upperbounds
    opt.ftol_rel = 1e-4
    opt.maxtime = 1
    opt.min_objective = (x,g) -> fittingobjective(x, g, model, syntheticdata, DPpoints, noise, baseparams)
    (minf,minx,ret) = optimize(opt, initguess)
    return (minf,minx,ret)
end

function multistartoptimizeparams(model, syntheticdata, DPpoints, noise, baseparams; npoints = 1000)
    if typeof(model)<:InitiationLimitedModel
        nparams = 1
        initguess = [baseparams.k_i]
        upperbounds = [1e5]
    end
    if typeof(model)<:Union{InitiationElongationModel, BasicTASEPModel}
        nparams = 2
        initguess = [baseparams.k_i,α]
        upperbounds = [1e5,1e3]
    end
    if typeof(model)<:LPTASEPModel
        nparams = 3
        initguess = [baseparams.k_i,α,γ]
        upperbounds = [1e5,1e3,1e2]
    end
    lowerbounds = zeros(nparams)
    (minf,minx,ret) = localoptimizeparams(model, syntheticdata, DPpoints, noise, baseparams, initguess, nparams, upperbounds, lowerbounds)
    for i in 1:npoints
        initguess = [rand(Uniform(lowerbounds[i],upperbounds[i])) for i in 1:nparams]
        (newminf,newminx,newret) = localoptimizeparams(model, syntheticdata, DPpoints, noise, baseparams, initguess, nparams, upperbounds, lowerbounds)
        if !isnan(newminf) && newminf<minf
            (minf,minx,ret) = (newminf,newminx,newret)
        end
    end
    cov = getcovariancematrix(model,minx, DPpoints, baseparams, noise)
    println("got $minf at $minx (returned $ret)")
    return (generatefullparameters(model,minx,baseparams),minx,cov)
end

#Take list of DNA inputs and list of polymerase inputs, outputs a list of DNA,polymerase input pairs for a full grid experimental design
function generategridinputs(DNApoints,Ppoints)
    (hcat([hcat([[D,P] for P in Ppoints]...) for D in DNApoints]...)')
end

function plotsyndata_constantD!(plt,syntheticdata, DPinputs, noise; color = :greens, modelname="", normalizexaxis = true, normalizeyaxis = true)
    DNApoints = unique(DPinputs[:,1])
    colorpal = palette(color, length(DNApoints))
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
function plotmodelprediction!(plt,model,params, covariancematrix, p_base, DPinputs;
    npoints = 100, 
    color = :blues,
    modelname="",
    normalizexaxis = true, 
    normalizeyaxis = true, 
    maximumT7 = nothing, 
    extraDNA = [], 
    kwargs...)

    DNAs = vcat(unique(DPinputs[:,1]), extraDNA)
    colorpal = palette(color, length(DNAs))
    if normalizexaxis
        if isnothing(maximumT7)
            maxT7 = maximum(relativeT7input)
        else
            maxT7 = maximumT7
        end
        xlabel = "T7 RNA Polymerase per DNA"
        relativeT7input = DPinputs[:,2] ./DPinputs[:,1]
        prange = (1e-5,maxT7)
    else
        if isnothing(maximumT7)
            maxT7 = maximum(DPinputs[:,2])
        else
            maxT7 = maximumT7
        end
        xlabel = "T7 RNA Polymerase (nM)"
        prange = (1e-5,maxT7)
    end
    if normalizeyaxis
        ylabel = "Transcription Rate (mM NTP/ hr) per nM DNA"
    else
        ylabel = "Transcription Rate (mM NTP/ hr)"
    end
    plot!(plt, xlabel = xlabel, ylabel = ylabel)
    for (ind,D) in enumerate(DNAs)
        plotQSrate_constantD!(plt, model, params, covariancematrix, p_base, D, prange; modellabel = modelname*": DNA = "*string(round(D, digits = 1)), npoints = npoints, color = colorpal[ind],normalizexaxis = normalizexaxis,normalizeyaxis = normalizeyaxis, kwargs...)
    end
    return plt
end

function generateplotsyntheticdata(model,params,DPinputs; noise = 0, kwargs...)
    plt = plot()
    plotmodelprediction!(plt,model,params,DPinputs; color = :greens, kwargs...)
    syntheticdata = generatesyntheticdata(model,params,DPinputs,noise)
    plotsyndata_constantD!(plt,syntheticdata, DPinputs, noise; kwargs...)
    return (plt,syntheticdata,DPinputs)
end


function plotparametricellipseprojection(parameters,parametercovariance;α = 0.05,labels = ["",""])
    χ²ₚ = Distributions.Chisq(3) # chi-squared distribution parameterized by p d.f.
    crit_χ² = Distributions.quantile(χ²ₚ, 1-α)
    crit_l = sqrt(crit_χ²)
    λ_covb, e_vec_covb = eigen(inv(parametercovariance)) # eigendecomposition of covariance
    t = collect(0.0:0.01:2.0*π)
    parametric_circle = vcat(cos.(t)',sin.(t)')
    major_minor_axes = crit_l * e_vec_covb * Diagonal(abs.(λ_covb).^(-0.5)) # map unit ball (ψ) to 2D ellipse confidence region (cov(b))
    parametric_ellipse = parameters .+major_minor_axes * parametric_circle
    plt = plot(parametric_ellipse[1,:],parametric_ellipse[2,:],xlabel = labels[1],ylabel = labels[2],label = "")
    return plt
end

function plotparametricellipse(parameters,parametercovariance; title = "", α = 0.05,labels = [L"k_i (h)","α","γ"])
    if length(parameters) ==1
        relativeuncertainty = 100*sqrt(parametercovariance[1])/parameters[1]
        println("Initiation Limited Model: "*string(round(relativeuncertainty,digits = 2))*" percent relative uncertainty")
    elseif length(parameters) == 2
        plt = plotparametricellipseprojection(parameters,parametercovariance;α = α,labels = labels[1:2])
        plot!(title = title, size = (300,300))
        return plt
    elseif length(parameters) == 3
        plots = []
        for i in 1:3
            p_reduced = parameters[Not([i])]
            cov_reduced = parametercovariance[Not([i]),Not([i])]
            labels_reducted = labels[Not([i])]
            plot_reduced = plotparametricellipseprojection(p_reduced,cov_reduced;α = α,labels = labels_reducted)
            plots = vcat(plots,plot_reduced)
        end
        return plot(plots...,layout = (1,3), size = (900,300),bottommargin = 5mm,leftmargin = 5mm,plot_title = title)
    end
end