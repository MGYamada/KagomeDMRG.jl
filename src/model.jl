"""
    spin_sites(lattice)

Create spin-1/2 ITensor indices conserving integer `Sz = 2Sᶻ`.
Reuse these indices when comparing or continuing states at different fluxes.
"""
spin_sites(lattice::KagomeCylinder) = siteinds("S=1/2", nsites(lattice); conserve_qns=true)

function _check_sites(sites, N)
    length(sites) == N || throw(ArgumentError("site count does not match the lattice"))
    all(s -> hasqns(s) && dim(s) == 2 && hastags(s, "S=1/2"), sites) ||
        throw(ArgumentError("expected U(1)-conserving spin-1/2 site indices"))
    for s in sites
        qn(s => 1) == QN("Sz", 1) && qn(s => 2) == QN("Sz", -1) ||
            throw(ArgumentError("site indices must use the Sz = 2Sᶻ convention"))
    end
    return nothing
end

"""
    twisted_exchange_mpo(sites, lattice, theta; gauge=:seam, hz=nothing)

Build a `ComplexF64` XXZ MPO. Each oriented bond contributes
`Jz Sz_i Sz_j + Jxy/2 (exp(im*A) S+_i S-_j + exp(-im*A) S-_i S+_j)`.
In seam gauge `A = wy*theta`; the longitudinal term is never twisted.
The optional vector `hz` adds `-sum(hz[i]*Sz_i)`. Nonuniform fields change
the model and are intended for explicitly labeled controls.
"""
function twisted_exchange_mpo(sites, lattice::KagomeCylinder, theta::Real;
                              gauge::Symbol=:seam, hz=nothing)
    N = nsites(lattice)
    _check_sites(sites, N)
    isfinite(theta) || throw(ArgumentError("theta must be finite"))
    gauge in (:seam, :uniform) || throw(ArgumentError("unknown gauge: $gauge"))
    fields = hz === nothing ? zeros(N) : Float64.(hz)
    length(fields) == N || throw(ArgumentError("hz must have one entry per site"))
    all(isfinite, fields) || throw(ArgumentError("hz must be finite"))
    terms = OpSum()
    nonzero = false
    for b in lattice.bonds
        if !iszero(b.Jz)
            terms += b.Jz, "Sz", b.i, "Sz", b.j
            nonzero = true
        end
        if !iszero(b.Jxy)
            a = (b.Jxy / 2) * cis(bond_phase(lattice, b, theta; gauge))
            terms += a, "S+", b.i, "S-", b.j
            terms += conj(a), "S-", b.i, "S+", b.j
            nonzero = true
        end
    end
    for i in eachindex(fields)
        if !iszero(fields[i])
            terms += -fields[i], "Sz", i
            nonzero = true
        end
    end
    if !nonzero
        # An empty OpSum has no well-defined MPO representation.
        return 0.0 * MPO(ComplexF64, sites, "Id")
    end
    return MPO(ComplexF64, terms, sites)
end
