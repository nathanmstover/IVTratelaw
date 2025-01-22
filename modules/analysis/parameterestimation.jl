function multistartoptimizeparams(model, syntheticdata, DPpoints, noise, baseparams; npoints = 1000, priorfunc = x->0, priormat = 0)
    α = 1+baseparams.N_all* baseparams.k_i/baseparams.k_e
    K = (baseparams.k_off+baseparams.k_i)/baseparams.k_on
    γ = baseparams.k_i*baseparams.f*baseparams.τ^2/(1+baseparams.f*baseparams.τ)
    if typeof(model)<:InitiationLimitedModel
        nparams = 2
        initguess = [baseparams.k_i,K]
        upperbounds = [1e5,1e3]
    end
    if typeof(model)<:InitiationElongationModel
        nparams = 3
        initguess = [baseparams.k_i,K,α]
        upperbounds = [1e5,1e3,1e3]
    end
    if typeof(model)<:BasicTASEPModel
        nparams = 4
        initguess = [baseparams.k_i,K,α,β]
        upperbounds = [1e5,1e3,1e3,1e3]
    end
    if typeof(model)<:LPTASEPModel
        nparams = 5
        initguess = [baseparams.k_i,K,α,β,γ]
        upperbounds = [1e5,1e3,1e3,1e3,1e3]
    end
    lowerbounds = zeros(nparams)
    (minf,minx,ret) = localoptimizeparams(model, syntheticdata, DPpoints, noise, baseparams, initguess, nparams, upperbounds, lowerbounds, priorfunc = priorfunc)
    for i in 1:npoints
        initguess = [rand(Uniform(lowerbounds[i],upperbounds[i])) for i in 1:nparams]
        (newminf,newminx,newret) = localoptimizeparams(model, syntheticdata, DPpoints, noise, baseparams, initguess, nparams, upperbounds, lowerbounds, priorfunc = priorfunc)
        if !isnan(newminf) && newminf<minf
            (minf,minx,ret) = (newminf,newminx,newret)
        end
    end
    cov = getcovariancematrix(model,minx, DPpoints, baseparams, noise, priormatrix = priormat)
    println("got $minf at $minx (returned $ret)")
    return (generatefullparameters(model,minx,baseparams),minx,cov)
end

function localoptimizeparams(model, syntheticdata, DPpoints, noise, baseparams, initguess, nparams, upperbounds, lowerbounds; kwargs...)
    opt = Opt(:LD_SLSQP, nparams)
    opt.lower_bounds = lowerbounds
    opt.upper_bounds = upperbounds
    opt.ftol_rel = 1e-4
    opt.maxtime = 1
    opt.min_objective = (x,g) -> fittingobjective(x, g, model, syntheticdata, DPpoints, noise, baseparams; kwargs...)
    (minf,minx,ret) = optimize(opt, initguess)
    return (minf,minx,ret)
end

function fittingobjective(x::Vector, grad::Vector, model, syntheticdata, DPpoints, noise, baseparams; priorfunc = x->0)
    if length(grad) > 0
        fun = p->modelresidual(model, p, syntheticdata, DPpoints, noise, baseparams; priorfunc = priorfunc)
        gnew = ForwardDiff.gradient(fun,x)
        for i in 1:length(grad)
            grad[i] = gnew[i]
        end
    end
    residual = modelresidual(model, x, syntheticdata, DPpoints, noise, baseparams; priorfunc = priorfunc)
    return residual
end

function modelresidual(model,estimatedparams, syntheticdata, DPpoints, noise, baseparams; priorfunc = x->0)
    predicteddata = modelprediction(model,estimatedparams, DPpoints, baseparams)
    priorfunc(estimatedparams)+norm((predicteddata .- syntheticdata) ./noise)^2#Sum of weighted squares equivalent to negative log likelihood
end







