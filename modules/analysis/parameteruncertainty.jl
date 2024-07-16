function getcovariancematrix(model,estimatedparams, DPpoints, baseparams, noise; priormatrix = 0)
    M = getinformationmatrix(model,estimatedparams, DPpoints, baseparams, noise, priormatrix = priormatrix)
    return inv(M)
end

function getinformationmatrix(model,estimatedparams, DPpoints, baseparams, noise; 
    priormatrix = zeros(length(estimatedparams),
    length(estimatedparams)))
    
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
    noise.^2 .* I(ndatapoints)
end
