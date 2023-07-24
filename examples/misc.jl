"""
function get_fe_exact(
    ::Val{:XYC},
    β::Float64,
    L::Int64;
    J::Float64 = 1.0
) -> free energy

    calculate the free energy for the isometric XY chain:
    H/J ~ SxSx + SySy

    ref:
    PRX 8, 031082 (2018)
    Eq. F3 and F4
"""
function get_fe_exact(
    β::Float64,
    L::Int64
)
    ε(k) = cos(k*π/(L+1))
    return -1/β * sum([log(1+exp(-β*ε(k))) for k in 1:L])
end
