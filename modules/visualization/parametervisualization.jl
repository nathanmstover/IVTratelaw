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

function plotparametricellipse(parameters,parametercovariance; title = "", α = 0.05,labels = [L"k_i (h)",L"K_{MD} (nM)","α","β","γ"])
    nparams = length(parameters)
    if nparams == 1
        relativeuncertainty = 100*sqrt(parametercovariance[1])/parameters[1]
        println("Initiation Limited Model: "*string(round(relativeuncertainty,digits = 2))*" percent relative uncertainty")
    elseif nparams > 1
        plots = []
        plotpairs = multiset_combinations(1:nparams,2)
        for pair in plotpairs
            p_reduced = parameters[pair]
            cov_reduced = parametercovariance[pair,pair]
            labels_reducted = labels[pair]
            plot_reduced = plotparametricellipseprojection(p_reduced,cov_reduced;α = α,labels = labels_reducted)
            plots = vcat(plots,plot_reduced)
        end
        nplots = length(plots)
        ncol = min(nplots,3)
        nrow = Int(ceil(nplots/3))
        return plot(plots...,layout = (nrow,ncol), size = (300*ncol,300*nrow), bottommargin = 5mm,leftmargin = 5mm,plot_title = title)
    end
end