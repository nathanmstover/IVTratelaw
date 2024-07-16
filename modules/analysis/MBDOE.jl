function multistartMBDOEoptimization(model, params, baseparams, nexperiments;
    npoints = 1000,
    kwargs...)
    (minf,minx,ret) = [Inf,0,0]
    for i in 1:npoints
        iterguess = reshape(choosefeasiblepoint(nexperiments;kwargs...),2*nexperiments)
        (newminf,newminx,newret) = localMBDOEoptimization(model, params, iterguess, baseparams; kwargs...)
        if !isnan(newminf) && newminf<minf
            (minf,minx,ret) = (newminf,newminx,newret)
        end
    end
    println("got $minf at $minx (returned $ret)")
    return (minf,minx,ret)
end

function choosefeasiblepoint(nexperiments;    
    upperDNA = 28.5, #nM
    lowerDNA = 1.0, #nM
    upperP = 384.0, #nM
    lowerP = 40.0, #nM
    DNAconcentration = 504,#stock solution nM
    T7concentration = 3200,#stock solution nM
    maxvolumeratio = 0.28,
    priormatrix = 0)#fraction of vial)

    pointls = []
    pointsaddedcounter = 0
    while pointsaddedcounter<nexperiments
        DNAadded = rand(Uniform(lowerDNA,upperDNA))
        T7added = rand(Uniform(lowerP,upperP))
        if DNAadded/DNAconcentration + T7added/T7concentration < maxvolumeratio
            append!(pointls,[[DNAadded,T7added]])
            pointsaddedcounter+=1
        end
    end
    return hcat(pointls...)'
end

function localMBDOEoptimization(model, params, initDPpoints, baseparams; 
    noise = (T7,DNA,param,base_params)-> 1.0, 
    upperDNA = 28.5, #nM
    lowerDNA = 1.0, #nM
    upperP = 384.0, #nM
    lowerP = 40.0, #nM
    DNAconcentration = 504,#stock solution nM
    T7concentration = 3200,#stock solution nM
    maxvolumeratio = 0.28,#fraction of vial
    priormatrix = zeros(length(params),length(params)))

    nvariables = Int(length(initDPpoints)/2)
    opt = Opt(:LD_SLSQP, nvariables*2)
    opt.lower_bounds = vcat(ones(nvariables)*lowerDNA,ones(nvariables)*lowerP)
    opt.upper_bounds = vcat(ones(nvariables)*upperDNA,ones(nvariables)*upperP)
    opt.ftol_rel = 1e-4
    opt.maxtime = 1#seconds
    inequality_constraint!(opt, (res,x,g) -> maxvolumeconstraint(res,x,g,DNAconcentration,T7concentration,maxvolumeratio), 1e-3*ones(nvariables))
    opt.min_objective = (x,g) -> Doptimalitywrapper(model, params, x, baseparams, noise, g, priormatrix = priormatrix)
    (minf,minx,ret) = optimize(opt, initDPpoints)
    return (minf,minx,ret)
end

function maxvolumeconstraint(result::Vector, x::Vector, grad::Matrix, DNAconc, T7conc, maxvolumeratio)#DNAconc and T7conc in nM, maxvolume in uL
    if length(grad) > 0
        for i in 1:Int(length(x)/2)
            grad[1,i] = 1/DNAconc
            grad[2,i] = 1/T7conc
        end
    end
    for i in 1:Int(length(x)/2)
        DNAadded = x[i]#nM
        T7added = x[i+Int(length(x)/2)]#nM
        result[i] = DNAadded/DNAconc + T7added/T7conc - maxvolumeratio
    end
end

function Doptimalitywrapper(model,estimatedparams, DPpointsflattened, baseparams, noise, grad; kwargs...)
    if length(grad) > 0
        fun = x->getDoptimalitymeasure(model,estimatedparams, x, baseparams, noise; kwargs...)
        gnew = ForwardDiff.gradient(fun,DPpointsflattened)
        for i in 1:length(grad)
            grad[i] = gnew[i]
        end
    end
    dmeasure = getDoptimalitymeasure(model,estimatedparams, DPpointsflattened, baseparams, noise; kwargs...)
    return dmeasure
end

function getDoptimalitymeasure(model,estimatedparams, DPpointsflattened, baseparams, noise; 
    priormatrix = zeros(length(estimatedparams),length(estimatedparams)))

    nexperiments = Int(length(DPpointsflattened)/2)
    DPpoints = reshape(DPpointsflattened,(nexperiments,2))
    noisevals = [noise(DPpoints[exptind,1],DPpoints[exptind,2],estimatedparams,baseparams) for exptind in 1:nexperiments]
    M = getinformationmatrix(model,estimatedparams, DPpoints, baseparams, noisevals, priormatrix = priormatrix)
    return det(inv(M))
end

