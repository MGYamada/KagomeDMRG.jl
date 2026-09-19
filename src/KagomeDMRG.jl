module KagomeDMRG

using LinearAlgebra
using KrylovKit
using OMEinsum
using Random
using SHA
using Serialization
using TOML
using ITensors
using ITensorMPS

include("provenance.jl")
include("lattice.jl")
include("extended_model.jl")
include("model.jl")
include("dmrg.jl")
include("observables.jl")
include("chirality.jl")
include("schmidt.jl")
include("checkpoint.jl")
include("continuation.jl")

export Bond, KagomeSite, KagomeCylinder, kagome_cylinder, nsites, site_index
export kagome_j1j2j3_cylinder, bond_families
export target_sector, right_region, reverse_bond, bond_phase, gauge_angles
export spin_sites, initial_mps, twisted_exchange_mpo, run_dmrg
export sz_profile, spin_correlations, spin_transfer, bond_energies
export KagomeTriangle, oriented_triangles, scalar_chirality_mpo, triangle_chiralities
export schmidt_diagnostics
export save_checkpoint, load_checkpoint, resume_dmrg
export FluxPolicy, continue_flux

function __init__()
    _initialize_execution_identity!()
    return nothing
end

end
