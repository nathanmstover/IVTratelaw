function getcovariancematrix(model,estimatedparams, DPpoints, baseparams, noise; priormatrix = 0)
    M = getinformationmatrix(model,estimatedparams, DPpoints, baseparams, noise, priormatrix = priormatrix)
    return inv(M)
end

function getinformationmatrix(model,estimatedparams, DPpoints, baseparams, noise; 
    priormatrix = zeros(length(estimatedparams),length(estimatedparams)))
    V = getexperimentalvariancematrix(model,estimatedparams, DPpoints, baseparams, noise)
    S = getsensitivitymatrix(model,estimatedparams, DPpoints, baseparams)
    return S*inv(V)*S' .+ priormatrix
end

function getsensitivitymatrix(model,estimatedparams, DPpoints, baseparams)
    fun = p->modelprediction(model,p, DPpoints, baseparams)
    ForwardDiff.jacobian(fun,estimatedparams)'
end

function getexperimentalvariancematrix(model,estimatedparams, DPpoints, baseparams, noise)
    ndatapoints = size(DPpoints)[1]
    (noise.^2/3) .* I(ndatapoints)#We need to convert experimental noise to variance of point. We do this by dividing by number of replicates.
end

function printuncertaintyestimates(p_thin,cov,N_all;nmc = 100000)
    nparams = length(p_thin)
    paramnames = ["k_i", "K_MD", "α", "β", "γ"] 
    ke_names = ["k_e,eff", "k_e,bp"]
    mcensemble = zeros(nmc,nparams)
    ke_ensemble = zeros(nmc,2)
    mean = p_thin
    if nparams>2
        ke_mean = [mean[1]/(mean[3]-1),N_all*mean[1]/(mean[3]-1)]
    end
    d = MvNormal(mean, Hermitian(cov))
    for i in 1:nmc
        x = rand(d, 1)
        mcensemble[i,:] = x
        if nparams>2
            ke_ensemble[i,1] = x[1]/(x[3]-1)
            ke_ensemble[i,2] = N_all*x[1]/(x[3]-1)
        end
    end
    lower_pointwise_CB = [percentile(mcensemble[:,j],100*0.05/2) for j in 1:nparams]
    upper_pointwise_CB = [percentile(mcensemble[:,j],100*(1-0.05/2)) for j in 1:nparams]
    lower_pointwise_CB_ke = [percentile(ke_ensemble[:,j],100*0.05/2) for j in 1:2]
    upper_pointwise_CB_ke = [percentile(ke_ensemble[:,j],100*(1-0.05/2)) for j in 1:2]

    [print(paramnames[i]*"    ") for i in 1:nparams]
    if nparams>2 
        [print(name*"  ") for name in ke_names]
    end
    println("")

    [print(string(round(meanval,digits = 1))*"  ") for meanval in mean]
    if nparams>2
        [print(string(round(ke_meanval,digits = 1))*"  ") for ke_meanval in ke_mean]
    end
    println("")

    [print(string(round(lowerCB,digits = 1))*"  ") for lowerCB in lower_pointwise_CB]
    if nparams>2
        [print(string(round(ke_lowerCB,digits = 1))*"  ") for ke_lowerCB in lower_pointwise_CB_ke]
    end
    println("")    
    
    [print(string(round(upperCB,digits = 1))*"  ") for upperCB in upper_pointwise_CB]
    if nparams>2
        [print(string(round(ke_upperCB,digits = 1))*"  ") for ke_upperCB in upper_pointwise_CB_ke]
    end
    println("")    
end
