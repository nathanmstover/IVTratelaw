function modelprediction(model,estimatedparams, DPpoints, baseparams)
    fullparams = generatefullparameters(model,estimatedparams,baseparams)
    generatesyntheticdata(model,fullparams,DPpoints,zeros(size(DPpoints)[1]))
end

function generatesyntheticdata(model,params,DPinputs,noise)
    ndatapoints = size(DPinputs)[1]
    syntheticdata = []
    for ind in 1:ndatapoints
        (DNA,T7) = DPinputs[ind,:]
        syntheticpoint = quasisteadyrate_analytic(model, params, DNA, T7) + noise[ind]*randn()
        append!(syntheticdata,syntheticpoint)
    end
    return syntheticdata
end

#Take list of DNA inputs and list of polymerase inputs, outputs a list of DNA,polymerase input pairs for a full grid experimental design
function generategridinputs(DNApoints,Ppoints)
    (hcat([hcat([[D,P] for P in Ppoints]...) for D in DNApoints]...)')
end



