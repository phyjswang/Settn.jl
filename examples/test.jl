include("../src/Settn.jl")
using .Settn
using .Settn.ITensors
using UnicodePlots

include("misc.jl")

let
    # XY chian
    L = 32
    s = siteinds("S=1/2", L; conserve_qns = false)
    os = OpSum()
    for j = 1:(L-1)
        os += 0.5,"S+",j,"S-",j+1
        os += 0.5,"S-",j,"S+",j+1
    end
    H = MPO(os, s)

    lsβ = [2.0^i for i in -12:-1]
    lsfe = zeros(length(lsβ))
    lsfe0 = copy(lsfe)
    for (βi, β) in enumerate(lsβ)
        println(prod(fill("=",42)))
        println("Start calculating ρ for β = $β, βi/nβ = $βi/$(length(lsβ))\n\n\n\n")
        ρ = get_ρ(
            H,
            β/2;
            outputlevel = 2,
            stop_tol = 1e-9,
            maxdim = 128,
            cutoff = 1e-16
        )
        lsfe[βi] = get_fe(ρ, β)
        lsfe0[βi] = get_fe_exact(β, L)
    end

    lineplot(
        lsβ,
        abs.((lsfe-lsfe0)./lsfe0),
        xscale = :log2,
        yscale = :log10
    )
end
