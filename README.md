# SETTN.jl

A prototype of the series-expansion thermal tensor network (SETTN) algorithm [[PRB **95**, 161104](https://doi.org/10.1103/PhysRevB.95.161104)] is illustrated for the XY chain. Other models can be easily customized thanks to `ITensors.jl`'s `OpSum` function.

Note that this package is also an essential ingredient for more advanced thermal tensor network algorithms, such as [XTRG.jl](https:://github.com/phyjswang/XTRG.jl) and [tanTRG.jl](https://github.com/phyjswang/tanTRG.jl). In particular, **two-site variational MPO sum and product** discussed in [PRB **95**, 161104](https://doi.org/10.1103/PhysRevB.95.161104) and [PRX **8**, 031082](https://doi.org/10.1103/PhysRevX.8.031082) are implemented, which can be used in other modules via:

```julia
using SETTN: contract, +
```


## Installation

```julia
julia> ] add Settn
```

## Usage

To make it minimal, this package only exports two functions,

```julia
"""
function get_ρ(H::MPO, β::Float64; kwargs...) -> ρ::MPO

    calculate the density matrix ρ(β) using SETTN
"""
function get_ρ(H::MPO, β::Float64; kwargs...)
```

and

```julia
"""
function get_fe(ρ::MPO, β::Float64) -> fe::Float64

    calculate free energy at β given ρ(β/2) and β
"""
function get_fe(ρ::MPO, β::Float64)
```

### Example

Check the `examples/` directory for computing free energy of the XY chain:

```julia
julia> cd("examples")
julia> using Pkg
julia> Pkg.instantiate()
julia> include("test.jl")
```

## Advanced usage

Two-site variational MPO sum `+` and product `contract`, are extensions from [ITensors.jl](https://github.com/ITensor/ITensors.jl) with `::Algorithm"variational"`, more details can be found in `./src/mpo.jl`; also check the relevant part of source code in `ITensors.jl` package. Note the sweep scheme for variational sum and product follows the same one for `ITensors.jl`'s `dmrg` function.
## TODO

- Implement a more numerically stable MPO product, such as the one given in [PRB **102**, 035147](https://doi.org/10.1103/PhysRevB.102.035147)
