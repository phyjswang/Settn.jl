using ITensors

using ITensors:
    AbstractMPS,
    AbstractProjMPO,
    nsite,
    OneITensor,
    site_range,
    leftlim,
    rightlim,
    setleftlim!,
    setrightlim!,
    checkflux,
    @debug_check,
    @printf,
    orthocenter,
    @timeit_debug,
    timer,
    set_nsite!,
    # promote_itensor_eltype,
    # datatype,
    # adapt,
    replaceprime!

import ITensors:
    position!,
    replacebond!,
    copy,
    contract,
    +,
    set_nsite!,
    length

using ITensors.NDTensors:
    Algorithm,
    @Algorithm_str

"""
    This ProjMPO computes and stores the projection of an
    MPO `M` into a basis defined by an MPO `M1`, leaving a
    certain number of site indices of the MPO unprojected.
    Which sites are unprojected can be shifted by calling
    the `position!` method.

    Drawing of the network represented by a ProjMPO `P(H)`,
    showing the case of `nsite(P)==2` and `position!(P,M,4)`
    for an MPO `M`:

     -        -        -                          -        -        -
    | |      | |      | |        |        |      | |      | |      | |
    | o------+-o------+-o--------o--------o------+-o------+-o------+-o    M1
    | |      | |      | |        |        |      | |      | |      | |
    |  -     |  -     |  -                       |  -     |  -     |  -
    |   |    |   |    |   |                      |   |    |   |    |   |
    |   |    |   |    |   |                      |   |    |   |    |   |
     -  |     -  |     -  |                       -  |     -  |     -  |
      | |      | |      | |                        | |      | |      | |
      o-+------o-+------o-+---                  ---o-+------o-+------o |  M*
      | |      | |      | |                        | |      | |      | |
       -        -        -                          -        -        -
"""
function _makeL!(P::AbstractProjMPO, M::MPO, k::Int)::Union{ITensor,Nothing}
    ll = P.lpos
    if ll ≥ k
        P.lpos = k
        return nothing
    end
    ll = max(ll, 0)
    L = lproj(P)
    while ll < k
        L = (L * dag(M)[ll + 1]) * P.H[ll + 1]
        P.LR[ll + 1] = L
        ll += 1
    end
    P.lpos = k
    return L
end

function makeL!(P::AbstractProjMPO, M::MPO, k::Int)
    _makeL!(P, M, k)
    return P
end

function _makeR!(P::AbstractProjMPO, M::MPO, k::Int)::Union{ITensor,Nothing}
    rl = P.rpos
    if rl ≤ k
        P.rpos = k
        return nothing
    end
    N = length(P.H)
    rl = min(rl, N + 1)
    R = rproj(P)
    while rl > k
        R = (R * dag(M)[rl - 1]) * P.H[rl - 1]
        P.LR[rl - 1] = R
        rl -= 1
    end
    P.rpos = k
    return R
end

function makeR!(P::AbstractProjMPO, M::MPO, k::Int)
    _makeR!(P, M, k)
    return P
end

"""
    position!(P::ProjMPO, M::MPO, pos::Int)

    Given an MPO `M`, shift the projection of the
    MPO represented by the ProjMPO `P` such that
    the set of unprojected sites begins with site `pos`.
    This operation efficiently reuses previous projections
    of the MPO on sites that have already been projected.
    The MPO `M` must have compatible bond indices with
    the previous projected MPO tensors for this
    operation to succeed.
"""
function position!(P::AbstractProjMPO, M::MPO, pos::Int)
    makeL!(P, M, pos - 1)
    makeR!(P, M, pos + nsite(P))
    return P
end

function _varisum_sweeps(;
    nsweeps = 16,
    maxdim = typemax(Int),
    mindim = 1,
    cutoff = 1E-16,
    kwargs...
)
    sweeps = Sweeps(nsweeps)
    setmaxdim!(sweeps, maxdim...)
    setmindim!(sweeps, mindim...)
    setcutoff!(sweeps, cutoff...)
    setnoise!(sweeps, 0.0)
    return sweeps
end

"""
    NB! M20 is expected to be closer to M10+M20
"""
function +(
    M10::MPO, M20::MPO, sweeps::Sweeps;
    nsite = 2,
    SwpConvCrit = 1e-16,
    minswp = 5,
    kwargs...
)
    @debug_check begin
        # Debug level checks
        # Enable with ITensors.enable_debug_checks()
        checkflux(M10)
        checkflux(M20)
    end

    if !hassameinds(siteinds, M10, M20)
        error("In `+(::MPS/MPO...)`, the input `MPS` or `MPO` do not have the same site
        indices. For example, the site indices of the first site are $(siteinds.([M10, M20], 1))")
    end

    which_decomp::Union{String,Nothing} = get(kwargs, :which_decomp, nothing)
    svd_alg::String = get(kwargs, :svd_alg, "divide_and_conquer")
    outputlevel::Int = get(kwargs, :outputlevel, 0)

    M1 = copy(M10)
    if !isortho(M1) || orthocenter(M1) != 1
        M1 = orthogonalize!(M1, 1)
    end
    @assert isortho(M1) && orthocenter(M1) == 1
    M2 = copy(M20)
    setprime(linkinds, M2, 2)
    if !isortho(M2) || orthocenter(M2) != 2
        M2 = orthogonalize!(M2, 1)
    end
    @assert isortho(M2) && orthocenter(M2) == 1

    N = length(M1)
    M3 = deepcopy(M2)
    setprime!(linkinds,M3,3)

    PM1 = ProjMPO(M1)
    PM1 = set_nsite!(PM1, nsite)
    PM1 = position!(PM1, M3, 1)

    PM2 = ProjMPO(M2)
    PM2 = set_nsite!(PM2, nsite)
    PM2 = position!(PM2, M3, 1)

    checknrm = zeros(nsweep(sweeps))
    for sw in 1:nsweep(sweeps)
        sw_time = @elapsed begin
            maxtruncerr = 0.0
            for (b, ha) in sweepnext(N)
                @debug_check begin
                    checkflux(M1)
                    checkflux(PM1)
                    checkflux(M2)
                    checkflux(PM2)
                end

                @timeit_debug timer "VariSum: position!" begin
                    PM1 = position!(PM1, M3, b)
                    PM2 = position!(PM2, M3, b)
                end

                @debug_check begin
                    checkflux(M1)
                    checkflux(PM1)
                    checkflux(M2)
                    checkflux(PM2)
                end

                ortho = ha == 1 ? "left" : "right"
                CA = contract(PM1)
                CB = contract(PM2)
                CE = CA + CB
                spec = replacebond!(
                    M3,
                    b,
                    CE;
                    maxdim=maxdim(sweeps, sw),
                    mindim=mindim(sweeps, sw),
                    cutoff=cutoff(sweeps, sw),
                    eigen_perturbation=nothing,
                    ortho=ortho,
                    which_decomp=which_decomp,
                    svd_alg=svd_alg,
                )
                maxtruncerr = max(maxtruncerr, spec.truncerr)

                @debug_check begin
                    checkflux(M1)
                    checkflux(PM1)
                    checkflux(M2)
                    checkflux(PM2)
                end

                if outputlevel >= 4
                    @printf("Sweep %d, half %d, bond (%d,%d) \n", sw, ha, b, b + 1)
                    @printf(
                      "  Truncated using cutoff=%.1E maxdim=%d mindim=%d\n",
                      cutoff(sweeps, sw),
                      maxdim(sweeps, sw),
                      mindim(sweeps, sw)
                    )
                    @printf(
                      "  Trunc. err=%.2E, bond dimension %d\n", spec.truncerr, dim(linkind(M3, b))
                    )
                    flush(stdout)
                end
            end
            checknrm[sw] = norm(M3)
        end
        if outputlevel >= 3
            @printf(
              "After sweep %d, maxlinkdim=%d maxerr=%.2E time=%.3f\n",
              sw,
              maxlinkdim(M3),
              maxtruncerr,
              sw_time
            )
            flush(stdout)
        end
        isdone = false
        if sw > minswp
            isdone = (checknrm[sw] - checknrm[sw-1])/checknrm[sw-1] < SwpConvCrit
        end
        if isdone && outputlevel >= 2
            println(
                "Variational MPO sum converged at $sw, with maxlinkdim=$(maxlinkdim(M3)) and maxerr=$maxtruncerr"
            )
        end
        isdone && break
        if sw == nsweep(sweeps)
            @warn "Variational MPO sum does NOT converge!"
        end
    end
    return noprime(linkinds,M3)
end

function +(::Algorithm"variational", x1::MPO, x2::MPO; kwargs...)
    return +(x1, x2, _varisum_sweeps(; kwargs...); kwargs...)
end

function contract(P::AbstractProjMPO)
    itensor_map = Union{ITensor,OneITensor}[lproj(P)]
    append!(itensor_map, P.H[site_range(P)])
    push!(itensor_map, rproj(P))

    if dim(first(itensor_map)) == 1
        reverse!(itensor_map)
    end

    Hv = ITensor(1.0)
    for it in itensor_map
        Hv *= it
    end
    return Hv
end

function replacebond!(M::MPO, b::Int, phi::ITensor; kwargs...)
    ortho::String = get(kwargs, :ortho, "left")
    swapsites::Bool = get(kwargs, :swapsites, false)
    which_decomp::Union{String,Nothing} = get(kwargs, :which_decomp, nothing)
    normalize::Bool = get(kwargs, :normalize, false)

    indsMb = inds(M[b])
    if swapsites
        sb = siteinds(M, b)
        sbp1 = siteinds(M, b+1)
        indsMb = replaceinds(indsMb, sb, sbp1)
    end

    L, R, spec = ITensors.factorize(
            phi,
            indsMb;
            which_decomp=which_decomp,
            tags=tags(linkind(M, b)),
            kwargs...
    )

    M[b] = L
    M[b+1] = R
    if ortho == "left"
        leftlim(M) == b - 1 && setleftlim!(M, leftlim(M) + 1)
        rightlim(M) == b + 1 && setrightlim!(M, rightlim(M) + 1)
        normalize && (M[b + 1] ./= norm(M[b + 1]))
    elseif ortho == "right"
        leftlim(M) == b && setleftlim!(M, leftlim(M) - 1)
        rightlim(M) == b + 2 && setrightlim!(M, rightlim(M) - 1)
        normalize && (M[b] ./= norm(M[b]))
    else
        error(
          "In replacebond!, got ortho = $ortho, only currently supports `left` and `right`."
        )
    end
    return spec
end

"""
    This Proj2MPOs computes and stores the projection of an
    MPO `M` into a basis defined by the product of two MPOs `Ma` and `Mb`,
    leaving a certain number of site indices of the MPO unprojected.
    Which sites are unprojected can be shifted by calling
    the `position!` method.

    Drawing of the network represented by a Proj2MPOs `P(Ma, Mb)`,
    showing the case of `nsite(P)==2` and `position!(P,M,4)`
    for an MPO `M`:

     -        -        -                          -        -        -
    | |      | |      | |        |        |      | |      | |      | |
    | o------+-o------+-o--------o--------o------+-o------+-o------+-o    Ma'
    | |      | |      | |        |        |      | |      | |      | |
    | o------+-o------+-o--------o--------o------+-o------+-o------+-o    Mb
    | |      | |      | |        |        |      | |      | |      | |
    |  -     |  -     |  -                       |  -     |  -     |  -
    |   |    |   |    |   |                      |   |    |   |    |   |
    |   |    |   |    |   |                      |   |    |   |    |   |
    |   |    |   |    |   |                      |   |    |   |    |   |
     -  |     -  |     -  |                       -  |     -  |     -  |
      | |      | |      | |                        | |      | |      | |
      o-+------o-+------o-+---                  ---o-+------o-+------o |  M*
      | |      | |      | |                        | |      | |      | |
       -        -        -                          -        -        -
"""
mutable struct Proj2MPOs <: AbstractProjMPO
    lpos::Int
    rpos::Int
    nsite::Int
    Ma::MPO
    Mb::MPO
    LR::Vector{ITensor}
end

Proj2MPOs(Ma::MPO, Mb::MPO, nsite::Int) = Proj2MPOs(0, length(Ma) + 1, nsite, Ma, Mb, Vector{ITensor}(undef, length(Ma)))

copy(P::Proj2MPOs) = Proj2MPOs(P.lpos, P.rpos, P.nsite, copy(P.Ma), copy(P.Mb), copy(P.LR))

length(P::Proj2MPOs) = length(P.Ma)

function set_nsite!(P::Proj2MPOs, nsite)
    P.nsite = nsite
    return P
end

function _makeL!(P::Proj2MPOs, M::MPO, k::Int)::Union{ITensor,Nothing}
    ll = P.lpos
    if ll ≥ k
        P.lpos = k
        return nothing
    end
    ll = max(ll, 0)
    L = lproj(P)
    while ll < k
        L = ((L * P.Ma[ll + 1]) * P.Mb[ll + 1]) * dag(M)[ll + 1]
        P.LR[ll + 1] = L
        ll += 1
    end
    P.lpos = k
    return L
end

function makeL!(P::Proj2MPOs, M::MPO, k::Int)
    _makeL!(P, M, k)
    return P
end

function _makeR!(P::Proj2MPOs, M::MPO, k::Int)::Union{ITensor,Nothing}
    rl = P.rpos
    if rl ≤ k
        P.rpos = k
        return nothing
    end
    N = length(P.Ma)
    rl = min(rl, N + 1)
    R = rproj(P)
    while rl > k
        R = ((R * P.Ma[rl - 1]) * P.Mb[rl - 1]) * dag(M)[rl - 1]
        P.LR[rl - 1] = R
        rl -= 1
    end
    P.rpos = k
    return R
end

function makeR!(P::Proj2MPOs, M::MPO, k::Int)
    _makeR!(P, M, k)
    return P
end

function contract(P::Proj2MPOs)::ITensor
    itensor_map = Union{ITensor,OneITensor}[lproj(P)]
    for i in site_range(P)
        push!(itensor_map, P.Ma[i])
        push!(itensor_map, P.Mb[i])
    end
    push!(itensor_map, rproj(P))

    if dim(first(itensor_map)) == 1
        reverse!(itensor_map)
    end

    Hv = ITensor(1.0)
    for it in itensor_map
        Hv *= it
    end
    return Hv
end

"""
    NB! M20 is expected to be closer to M10*M20
"""
function contract(M10::MPO, M20::MPO, sweeps::Sweeps;
    nsite = 2,
    SwpConvCrit = 1e-16,
    minswp = 5,
    kwargs...
)
    @debug_check begin
        # Debug level checks
        # Enable with ITensors.enable_debug_checks()
        checkflux(M10)
        checkflux(M20)
    end

    if !hassameinds(siteinds, replaceprime(M10, 2=>0), M20)
        error("In `+(::MPS/MPO...)`, the input `MPS` or `MPO` do not have the same site
        indices. For example, the site indices of the first site are $(siteinds.([M10, M20], 1))")
    end

    which_decomp::Union{String,Nothing} = get(kwargs, :which_decomp, nothing)
    svd_alg::String = get(kwargs, :svd_alg, "divide_and_conquer")
    outputlevel::Int = get(kwargs, :outputlevel, 0)

    M1 = copy(M10)
    if !isortho(M1) || orthocenter(M1) != 1
        M1 = orthogonalize!(M1, 1)
    end
    @assert isortho(M1) && orthocenter(M1) == 1
    M2 = copy(M20)
    setprime!(linkinds,M2,2)
    if !isortho(M2) || orthocenter(M2) != 1
        M2 = orthogonalize!(M2, 1)
    end
    @assert isortho(M2) && orthocenter(M2) == 1

    N = length(M1)
    M3 = deepcopy(M2)
    replaceprime!(M3, 1=>2; tags="Site")
    setprime!(linkinds, M3, 3)
    PM = Proj2MPOs(M1, M2, nsite)
    PM = position!(PM, M3, 1)
    checknrm = zeros(nsweep(sweeps))
    for sw in 1:nsweep(sweeps)
        sw_time = @elapsed begin
            maxtruncerr = 0.0

            for (b, ha) in sweepnext(N)
                @debug_check begin
                    checkflux(M1)
                    checkflux(M2)
                    checkflux(PM)
                end

                @timeit_debug timer "VariSum: position!" begin
                    PM = position!(PM, M3, b)
                end

                @debug_check begin
                    checkflux(M1)
                    checkflux(M2)
                    checkflux(PM)
                end

                ortho = ha == 1 ? "left" : "right"
                CE = contract(PM)
                spec = replacebond!(
                    M3,
                    b,
                    CE;
                    maxdim=maxdim(sweeps, sw),
                    mindim=mindim(sweeps, sw),
                    cutoff=cutoff(sweeps, sw),
                    eigen_perturbation=nothing,
                    ortho=ortho,
                    which_decomp=which_decomp,
                    svd_alg=svd_alg,
                )
                maxtruncerr = max(maxtruncerr, spec.truncerr)

                @debug_check begin
                    checkflux(M1)
                    checkflux(M2)
                    checkflux(M3)
                    checkflux(PM)
                end

                if outputlevel >= 4
                    @printf("Sweep %d, half %d, bond (%d,%d) \n", sw, ha, b, b + 1)
                    @printf(
                      "  Truncated using cutoff=%.1E maxdim=%d mindim=%d\n",
                      cutoff(sweeps, sw),
                      maxdim(sweeps, sw),
                      mindim(sweeps, sw)
                    )
                    @printf(
                      "  Trunc. err=%.2E, bond dimension %d\n", spec.truncerr, dim(linkind(M3, b))
                    )
                    flush(stdout)
                end
            end # for (b, ha) in sweepnext(N)
            checknrm[sw] = norm(M3)
        end # for sw in 1:nsweep(sweeps)
        if outputlevel >= 3
            @printf(
              "After sweep %d, maxlinkdim=%d maxerr=%.2E time=%.3f\n",
              sw,
              maxlinkdim(M3),
              maxtruncerr,
              sw_time
            )
            flush(stdout)
        end
        isdone = false
        if sw > minswp
            isdone = (checknrm[sw] - checknrm[sw-1])/checknrm[sw-1] < SwpConvCrit
        end
        if isdone && outputlevel >= 2
            println(
                "Variational MPO product converged at $sw, with maxlinkdim=$(maxlinkdim(M3)) and maxerr=$(maxtruncerr)"
            )
        end
        isdone && break
        if sw == nsweep(sweeps)
            @warn "Variational MPO sum does NOT converge!"
        end
    end
    return noprime(linkinds, M3)
end

function contract(::Algorithm"variational", x1::MPO, x2::MPO; kwargs...)
    return contract(x1, x2, _varisum_sweeps(; kwargs...); kwargs...)
end
