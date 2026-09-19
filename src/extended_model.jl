# Representatives on the infinite lattice. Each tuple is
# (family, source sublattice, destination sublattice, delta x, delta y).
# J3 connects opposite hexagon vertices, excluding straight-chain neighbors
# at the same distance. See docs/research/p3_csl_control_design.md.
const _KAGOME_EXTENDED_BOND_TEMPLATES = (
    (:J2, :A, :C, 1, -1),
    (:J2, :B, :A, 1, -1),
    (:J2, :B, :C, 1, 0),
    (:J2, :C, :A, 1, 0),
    (:J2, :A, :B, 0, -1),
    (:J2, :B, :C, 0, -1),
    (:J3, :A, :A, 1, -1),
    (:J3, :B, :B, 0, -1),
    (:J3, :C, :C, 1, 0),
)

function _foreach_extended_bond(visit, lattice)
    for x in 0:(lattice.Lx - 1), y in 0:(lattice.Ly - 1)
        for (family, source, destination, dx, dy) in _KAGOME_EXTENDED_BOND_TEMPLATES
            # Keep a bond whenever both endpoints survive the open boundary;
            # its parent hexagon need not survive as a complete plaquette.
            xj, yj = x + dx, y + dy
            0 <= xj < lattice.Lx || continue
            wy, wrapped_y = fldmod(yj, lattice.Ly)
            i = site_index(lattice, x, y, source)
            j = site_index(lattice, xj, wrapped_y, destination)
            visit(family, i, j, wy)
        end
    end
    return nothing
end

"""
    kagome_j1j2j3_cylinder(Lx, Ly; J1=1.0, J2=0.5, J3=0.5)

Build an isotropic J1-J2-J3 kagome cylinder with the same sites, ordering,
open-axis termination and periodic wrap as [`kagome_cylinder`](@ref).
J2 skips one vertex of a hexagon; J3 connects only opposite hexagon vertices,
not the distinct straight-chain neighbors at the same distance.

Define the bonds on the infinite lattice, retain those with both endpoints
inside the open cylinder, and record the signed image winding before wrapping
the periodic coordinate. This retains valid edge bonds of incomplete hexagons.
The family counts are `(6Lx-2)Ly`, `(6Lx-4)Ly`, and `(3Lx-2)Ly` for J1, J2,
and J3. Zero-coupling bonds remain in the corresponding family explicitly.
All couplings must be finite real numbers; both `Jxy` and `Jz` equal `Ja`.
The nearest-neighbor builder and its defaults are unchanged.
"""
function kagome_j1j2j3_cylinder(Lx::Integer, Ly::Integer;
                               J1::Real=1.0, J2::Real=0.5, J3::Real=0.5)
    # Validate every coupling, including families with zero coefficient.
    couplings = (J1=Bond(1, 2, J1, J1, 0).Jxy,
                 J2=Bond(1, 2, J2, J2, 0).Jxy,
                 J3=Bond(1, 2, J3, J3, 0).Jxy)
    lattice = kagome_cylinder(Lx, Ly; Jxy=couplings.J1, Jz=couplings.J1)
    _foreach_extended_bond(lattice) do family, i, j, wy
        coupling = getproperty(couplings, family)
        push!(lattice.bonds, Bond(i, j, coupling, coupling, wy))
    end
    return lattice
end

_family_bond_key(i, j, wy) = i < j ? (i, j, wy) : (j, i, -wy)

"""
    bond_families(lattice) -> Vector{Symbol}

Classify each actual bond as `:J1`, `:J2`, `:J3`, or `:other` by its endpoints
and signed periodic image in the documented kagome geometry. Classification
is unchanged by reversing a bond or changing its couplings, including to zero.
J3 excludes straight-chain bonds. The lattice must use the canonical sites
and ordering of [`kagome_cylinder`](@ref); arbitrary additional bonds are allowed.

Build one geometric lookup in O(number of sites), then classify the supplied
bonds in O(number of bonds). No model identity or family labels are stored in
the lattice, so edits to its bond vector are reflected on every call.
"""
function bond_families(lattice::KagomeCylinder)::Vector{Symbol}
    reference = kagome_cylinder(lattice.Lx, lattice.Ly)
    lattice.sites == reference.sites ||
        throw(ArgumentError("bond families require the documented kagome sites and ordering"))
    lookup = Dict(_family_bond_key(b.i, b.j, b.wy) => :J1 for b in reference.bonds)
    _foreach_extended_bond(reference) do family, i, j, wy
        lookup[_family_bond_key(i, j, wy)] = family
    end
    return Symbol[get(lookup, _family_bond_key(b.i, b.j, b.wy), :other)
                  for b in lattice.bonds]
end
