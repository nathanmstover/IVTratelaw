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

function modelresidual(model,estimatedparams, syntheticdata, DPpoints, noise, baseparams)
    predicteddata = modelprediction(model,estimatedparams, DPpoints, baseparams)
    norm((predicteddata .- syntheticdata) ./noise)^2#Sum of weighted squares equivalent to negative log likelihood
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





