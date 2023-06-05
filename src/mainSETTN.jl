function getρ(H::MPO, s, β; kwargs...)
    nmax = get(kwargs, :nmax, 1024)
    exactmax = get(kwargs, :exactmax, 3)

    H0 = MPO(s, "Id")
    # nrm0 = tr(H0)
    nrm0 = norm(H0)
    H0 /= nrm0
    rho1 = H0

    Hn = copy(H)
    Hn /= nrm0

    rho1::MPO = +(rho1, -β * Hn; alg="directsum")

    lognrmtot = 0
    for i = 2:nmax
        if i < exactmax
            Hn::MPO = apply(H, Hn; alg="zipup", kwargs...)
        else
            Hn::MPO = apply(H, Hn; alg="variational", kwargs...)
        end
        nrm = norm(Hn)
        Hn /= nrm
        lognrmtot += log(nrm)
        coeff = Float64(((-1)^i*exp(lognrmtot + (i)*log(big(β)) - log(factorial(big(i))))))
        nrm1 = norm(rho1)
        nrm2 = coeff
        if i < exactmax
            rho1 = +(rho1, coeff*Hn; alg="directsum")
        else
            if nrm1 > nrm2
                rho1::MPO = +(coeff/nrm1*Hn, rho1/nrm1; alg="variational", kwargs...)
                rho1 *= nrm1
            else
                rho1::MPO = +(rho1/nrm2, Hn; alg="variational", kwargs...)
                rho1 *= nrm2
            end
        end
        nrmnew = norm(rho1)
        stopQ = abs((nrmnew-nrm1)/nrm1)
        if stopQ < 1e-16
            println("ρ converged at $i/$nmax")
            break
        end
        if i == nmax
            println("ρ NOT converged")
        end
    end
    return rho1, nrm0
end

function mainSETTN(H::MPO, s, lsβ, opnames::Vector{String}; kwargs...)
    lsex = []
    lsfe = Float64[]
    for (idx, β) in enumerate(lsβ)
        println("for β = $β, i = $idx / $(length(lsβ))")
        rho1, nrm0 = getρ(H, s, β/2; kwargs...)
        push!(lsfe, -1/β * (2*log(nrm0) + 2* log(norm(rho1))))
        push!(lsex,expect(rho1, opnames))
    end
    return lsfe, lsex
end
