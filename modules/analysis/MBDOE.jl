function multistartMBDOEoptimization(model, params, initDPpoints, baseparams;
    upperDNA = 28.5, 
    lowerDNA = 1.0, 
    upperP = 384.0, 
    lowerP = 40.0, 
    npoints = 1000,
    kwargs...)

    nvariables = Int(length(initDPpoints)/2)

    lowerbounds = vcat(ones(nvariables)*lowerDNA,ones(nvariables)*lowerP)
    upperbounds = vcat(ones(nvariables)*upperDNA,ones(nvariables)*upperP)
    (minf,minx,ret) = localMBDOEoptimization(model, params, initDPpoints, baseparams; upperDNA = upperDNA, lowerDNA = lowerDNA, upperP = upperP, lowerP = lowerP, kwargs...)

    for i in 1:npoints
        iterguess = [rand(Uniform(lowerbounds[j],upperbounds[j])) for j in 1:2*nvariables]
        (newminf,newminx,newret) = localMBDOEoptimization(model, params, iterguess, baseparams; upperDNA = upperDNA, lowerDNA = lowerDNA, upperP = upperP, lowerP = lowerP, kwargs...)
        if !isnan(newminf) && newminf<minf
            (minf,minx,ret) = (newminf,newminx,newret)
        end
    end
    println("got $minf at $minx (returned $ret)")
    return (minf,minx,ret)
end

function localMBDOEoptimization(model, params, initDPpoints, baseparams; 
    noise = (T7,DNA,param,base_params)-> 1.0, 
    upperDNA = 28.5, 
    lowerDNA = 1.0, 
    upperP = 384.0, 
    lowerP = 40.0, 
    priormatrix = zeros(length(params),length(params)))

    nvariables = Int(length(initDPpoints)/2)
    opt = Opt(:LD_SLSQP, nvariables*2)
    opt.lower_bounds = vcat(ones(nvariables)*lowerDNA,ones(nvariables)*lowerP)
    opt.upper_bounds = vcat(ones(nvariables)*upperDNA,ones(nvariables)*upperP)
    opt.ftol_rel = 1e-10
    opt.maxtime = 1#seconds
    opt.min_objective = (x,g) -> Doptimalitywrapper(model, params, x, baseparams, noise, g, priormatrix = priormatrix)
        (minf,minx,ret) = optimize(opt, initDPpoints)
    return (minf,minx,ret)
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
    priormatrix = zeros(length(estimatedparams),
    length(estimatedparams)))
    
    nexperiments = Int(length(DPpointsflattened)/2)
    DPpoints = reshape(DPpointsflattened,(nexperiments,2))
    noisevals = [noise(DPpoints[exptind,1],DPpoints[exptind,2],estimatedparams,baseparams) for exptind in 1:nexperiments]
    M = getinformationmatrix(model,estimatedparams, DPpoints, baseparams, noisevals, priormatrix = priormatrix)
    return det(inv(M))
end

