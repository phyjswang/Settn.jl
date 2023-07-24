using ITensors

"""
function get_ρ(H::MPO, β::Float64; kwargs...) -> ρ::MPO

    calculate the density matrix ρ(β) using SETTN
"""
function get_ρ(H::MPO, β::Float64; kwargs...)
    maxord = get(kwargs, :maxord, 1024)
    exactmax = get(kwargs, :exactmax, 3)
    stop_tol = get(kwargs, :stop_tol, 1e-12)
    outputlevel = get(kwargs, :outputlevel, 0)

    s = firstsiteinds(H)
    s = dag(s) # for symmetric Indices

    ρ = MPO(s, "Id")
    Hn = copy(H)
    ρ = +(ρ, -β * Hn; alg="directsum")
    lognrmρ1 = lognorm(ρ)
    for i = 2:maxord
        lognrmρ = lognrmρ1
        if i > exactmax
            Hn = apply(H, Hn; alg="variational", kwargs...)
        else
            Hn = apply(H, Hn; alg="zipup", kwargs...)
        end
        logcoef = i*log(β) - sum([log(i1) for i1 in 1:i])
        nrm1 = lognorm(ρ)
        nrm2 = logcoef + lognorm(Hn)
        if i > exactmax
            # NB! The second MPO is expected to be closer to the resulting MPO
            if nrm1 > nrm2
                ρ = +((-1)^i*exp(logcoef)*Hn/exp(nrm1), ρ/exp(nrm1); alg="variational", kwargs...)
                ρ *= exp(nrm1)
            else
                ρ = +(ρ/exp(nrm2), (-1)^i*exp(logcoef)*Hn/exp(nrm2); alg="variational", kwargs...)
                ρ *= exp(nrm2)
            end
        else
            ρ = +(ρ, (-1)^i*exp(logcoef)*Hn; alg="directsum")
        end
        lognrmρ1 = lognorm(ρ)
        reldiff = (lognrmρ-lognrmρ1)/lognrmρ
        if abs(reldiff) < stop_tol
            println("ρ converged at $i/$maxord")
            break
        else
            if outputlevel > 0
                println("Free energy rel. err. = $reldiff at n=$i")
                println(join(vcat(["-" for _ in 1:42],["\n"])))
            end
        end
        if i == maxord
            println("ρ NOT converged")
        end
    end
    return ρ
end

"""
function get_fe(ρ::MPO, β::Float64) -> fe::Float64

    calculate free energy at β given ρ(β/2) and β
"""
function get_fe(ρ::MPO, β::Float64)
    return -1/β * (2*lognorm(ρ))
end
