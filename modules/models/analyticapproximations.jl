#Takes parameters and returns general case for analytic IC
function quasisteadyrate_analytic_general(k_i, K, α, β, θ, γ, Dtot, polymerasetot)
    p = polymerasetot ./Dtot
    cubicterm = α .*γ
    quadraticterm  = α .+ β ./Dtot .-γ .*(p .+θ ./ Dtot .+α)
    linearterm = -(p .+ α .+K ./ Dtot .-γ .* p)
    constantterm = p
    ic = polynomialroots(constantterm, linearterm, quadraticterm, cubicterm)
    rate = Dtot .*k_i .*ic .*(1 .-(β .*ic)/(K .-θ))/(1 .+γ .*ic)#Rate in nM RNA strands per hour
    return rate, ic 
end

#Takes same arguments as runtranscriptionmodel, returns quasisteady rate of RNA production
function quasisteadyrate_analytic(model::TranscriptionModel, params::AbstractArray{T1}, DNA, P) where {T1<:Number}
    (α, β, θ, γ) = (1, 0, 0, 0)
    K = (params.k_i+params.k_off)/params.k_on
    if typeof(model)<:Union{InitiationElongationModel, ContinuumModel}
        α = params.α
    end
    if typeof(model)<:ContinuumModel
        β = (α-1)*params.k_i*params.L/(params.k_on*params.N_all)
        θ = params.k_off/params.k_on
    end
    if typeof(model)<:LPTASEPModel
        γ = params.γ
    end
    (rate, ic) = quasisteadyrate_analytic_general(params.k_i, K, α, β, θ, γ, DNA, P)

    if typeof(model)<:BasicTASEPModel && ((β*ic)/(K-θ) .>0.5)[1]#Over approximation threshold
        rate = params.k_i*DNA*(K-θ)/(4*β)
    end
    if typeof(model)<:LPTASEPModel
        if 1+γ*(K-θ)/β > 0
            amax = (sqrt(1+γ*(K-θ)/β)-1)/(γ*(K-θ)/β)
        else
            amax = 0
        end
        if ((β*ic)/(K-θ).>amax)[1]#Over approximation threshold
        rate = params.k_i*DNA*((K-θ)/β)*amax*(1-amax)/(1+amax*γ*((K-θ)/β))
        end
    end
    return rate*params.N_all*1e-6#Output in mM NTP per hour
end

function polynomialroots(d, c, b, a)
    #Use quadratic equation in case of no cubic term - faster?
    if a == 0.0
        return (-c .-sqrt.(c .^2 .-4 .*b .*d)) ./(2 .*b)
    else
        (a,b,c,d) = (Complex.(i) for i in [a,b,c,d])
        return solve_cubic_eq((d,c,b,a))
    end
end

function solve_cubic_eq(poly)#Adapted from PolynomialRoots.jl
    third = 1//3
    # Cubic equation solver for complex polynomial (degree=3)
    # http://en.wikipedia.org/wiki/Cubic_function   Lagrange's method
    a1  =  1 / poly[4]
    E1  = -poly[3] .*a1
    E2  =  poly[2] .*a1
    E3  = -poly[1] .*a1
    s0  =  E1
    E12 =  E1 .*E1
    A   =  2 .*E1 .*E12  .- 9 .*E1 .*E2  .+ 27 .*E3 # = s1^3 + s2^3
    B   =  E12 .- 3*E2                 # = s1 s2
    # quadratic equation: z^2 - Az + B^3=0  where roots are equal to s1^3 and s2^3
    Δ = sqrt.(A .*A .- 4 .*B .*B .*B)
    if (real.(conj.(A) .*Δ) .>=0)[1] # scalar product to decide the sign yielding bigger magnitude
        s1 = exp.(log.(0.5 .* (A .+ Δ)) .* third)
    else
        s1 = exp.(log.(0.5 .* (A .- Δ)) .* third)
    end
    if s1 == 0
        s2 = s1
    else
        s2 = B ./ s1
    end
    zeta1 = complex(-0.5, sqrt((3.0))*0.5)
    zeta2 = conj(zeta1)
    return real(third*(s0 + s1*zeta2 + s2*zeta1))
    # xs = [i for i in xstuple]
    # return real(xs[2])#Choosing physically relevant root
    # #return real((xs[0 .<real.(xs).<1 .&& abs.(imag.(xs)) .= 1e-10])[1])#Choosing physically relevant root
end
