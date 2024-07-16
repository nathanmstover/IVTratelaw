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