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