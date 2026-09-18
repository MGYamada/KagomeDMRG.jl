"""
    Bond(i, j, Jxy, Jz, wy)

One oriented exchange bond, counted once. The endpoint indices are positive
and distinct. `wy` is the signed periodic image of `j` relative to `i` along
the cylinder wrap vector. Reversing the orientation reverses `wy`.

In seam gauge the coefficient of `S+_i S-_j` is
`(Jxy / 2) * cis(wy * theta)`; the reverse exchange uses its conjugate.
"""
struct Bond
    i::Int
    j::Int
    Jxy::Float64
    Jz::Float64
    wy::Int

    function Bond(i::Integer, j::Integer, Jxy::Real, Jz::Real, wy::Integer)
        i > 0 && j > 0 || throw(ArgumentError("bond indices must be positive"))
        i != j || throw(ArgumentError("self bonds are not supported"))
        xy, z = Float64(Jxy), Float64(Jz)
        isfinite(xy) && isfinite(z) ||
            throw(ArgumentError("bond couplings must be finite Float64 values"))
        return new(Int(i), Int(j), xy, z, Int(wy))
    end
end

"""Return the same physical bond with reversed endpoints and winding."""
reverse_bond(b::Bond) = Bond(b.j, b.i, b.Jxy, b.Jz, -b.wy)

"""
A kagome site with a one-based MPS `index`, zero-based cell coordinates
`x, y`, sublattice `:A`, `:B`, or `:C`, and Cartesian `position`.
"""
struct KagomeSite
    index::Int
    x::Int
    y::Int
    sublattice::Symbol
    position::NTuple{2,Float64}
end

"""
A kagome cylinder open in cell coordinate `x` and periodic in `y`.
The primitive vectors are `a1=(1,0)` and `a2=(1/2,sqrt(3)/2)`;
the wrap vector is `Ly*a2`. Sites are ordered by `x`, then `y`, then `A,B,C`.
`Ly` counts primitive cells and is not a literature YC label.
"""
struct KagomeCylinder
    Lx::Int
    Ly::Int
    sites::Vector{KagomeSite}
    bonds::Vector{Bond}
end

"""Return the number of physical spin sites."""
nsites(lattice::KagomeCylinder) = length(lattice.sites)

function _sublattice_index(sublattice::Symbol)
    sublattice === :A && return 1
    sublattice === :B && return 2
    sublattice === :C && return 3
    throw(ArgumentError("sublattice must be :A, :B, or :C"))
end

"""
    site_index(lattice, x, y, sublattice)

Return the one-based site index of a canonical cell coordinate. Coordinates
must satisfy `0 <= x < Lx` and `0 <= y < Ly`; this function does not wrap them.
"""
function site_index(lattice::KagomeCylinder, x::Integer, y::Integer,
                    sublattice::Symbol)
    0 <= x < lattice.Lx || throw(ArgumentError("x is outside the open cylinder"))
    0 <= y < lattice.Ly || throw(ArgumentError("y must be a canonical periodic coordinate"))
    return 3 * (Int(x) * lattice.Ly + Int(y)) + _sublattice_index(sublattice)
end

"""
    kagome_cylinder(Lx, Ly; Jxy=1.0, Jz=1.0)

Build the nearest-neighbor kagome cylinder with `Lx >= 1` and `Ly >= 3`.
The sublattice offsets are `A=0`, `B=a1/2`, and `C=a2/2`. For each cell,
include the three intracell bonds and the bonds
`A(x,y)-B(x-1,y)`, `A(x,y)-C(x,y-1)`, and `B(x,y)-C(x+1,y-1)`.
Omit bonds whose open-coordinate endpoint lies outside the cylinder.

The result has `3Lx*Ly` sites and `(6Lx-2)*Ly` bonds. At either open end
only bonds leaving the retained cells are removed; no extra edge couplings
are added. Winding records the periodic image before wrapping `y`.
"""
function kagome_cylinder(Lx::Integer, Ly::Integer; Jxy::Real=1.0, Jz::Real=1.0)
    Lx >= 1 || throw(ArgumentError("Lx must be at least 1"))
    Ly >= 3 || throw(ArgumentError("Ly must be at least 3; narrower cylinders are not supported"))
    nx, ny = Int(Lx), Int(Ly)
    n = Base.checked_mul(3, Base.checked_mul(nx, ny))
    # Validate couplings even before constructing the first physical bond.
    coupling = Bond(1, 2, Jxy, Jz, 0)
    sites = KagomeSite[]
    sizehint!(sites, n)
    offsets = ((0.0, 0.0), (0.5, 0.0), (0.25, sqrt(3.0) / 4))
    labels = (:A, :B, :C)
    for x in 0:(nx - 1), y in 0:(ny - 1), sub in 1:3
        dx, dy = offsets[sub]
        position = (x + y / 2 + dx, y * sqrt(3.0) / 2 + dy)
        push!(sites, KagomeSite(length(sites) + 1, x, y, labels[sub], position))
    end
    bonds = Bond[]
    sizehint!(bonds, Base.checked_mul(6 * nx - 2, ny))
    lattice = KagomeCylinder(nx, ny, sites, bonds)
    function add_bond(xi, yi, si, xj, yj, sj)
        0 <= xj < nx || return nothing
        wy, wrapped_y = fldmod(yj, ny)
        i = site_index(lattice, xi, yi, si)
        j = site_index(lattice, xj, wrapped_y, sj)
        push!(bonds, Bond(i, j, coupling.Jxy, coupling.Jz, wy))
        return nothing
    end
    for x in 0:(nx - 1), y in 0:(ny - 1)
        add_bond(x, y, :A, x, y, :B)
        add_bond(x, y, :A, x, y, :C)
        add_bond(x, y, :B, x, y, :C)
        add_bond(x, y, :A, x - 1, y, :B)
        add_bond(x, y, :A, x, y - 1, :C)
        add_bond(x, y, :B, x + 1, y - 1, :C)
    end
    return lattice
end

"""
    target_sector(N)

Return `(N, Q, M, Nup, Ndown)` for `M/Msat = 1/9`. Here `Q=2M` is the
integer tensor charge and `M` is physical total `Sz`. Require a positive
number of sites divisible by nine; odd `N` correctly gives half-integer `M`.
"""
function target_sector(N::Integer)
    N > 0 || throw(ArgumentError("the number of sites must be positive"))
    N % 9 == 0 || throw(ArgumentError("the 1/9 magnetization target requires N divisible by 9"))
    n = Int(N)
    q = n ÷ 9
    return (N=n, Q=q, M=q / 2, Nup=5q, Ndown=4q)
end

"""
    right_region(lattice, cut)

Return the site indices in cells with `x >= cut`. The cut lies between
complete cell columns and must satisfy `1 <= cut < Lx`. Use the measured
zero-flux density on these same sites when defining spin transfer.
"""
function right_region(lattice::KagomeCylinder, cut::Integer)
    1 <= cut < lattice.Lx || throw(ArgumentError("cut must lie between interior cell columns"))
    return [site.index for site in lattice.sites if site.x >= cut]
end

function _finite_theta(theta::Real)
    value = Float64(theta)
    isfinite(value) || throw(ArgumentError("theta must be finite"))
    return value
end

"""
    gauge_angles(lattice, theta)

Return the local rotation angles `chi_i = -theta*eta_i/Ly`, where
`eta_i = y_i + (sublattice_i == :C ? 1/2 : 0)`. With
`U=exp(im*sum(chi_i*Sz_i))`, the uniform-gauge Hamiltonian is `U*Hseam*U'`.
"""
function gauge_angles(lattice::KagomeCylinder, theta::Real)
    value = _finite_theta(theta)
    return [-value * (site.y + (site.sublattice === :C ? 0.5 : 0.0)) / lattice.Ly
            for site in lattice.sites]
end

"""
    bond_phase(lattice, bond, theta; gauge=:seam)

Return the real Peierls angle `A_b`, not its exponential. The exchange
coefficient is `(bond.Jxy/2)*cis(A_b)`. In `:seam` gauge `A_b=wy*theta`;
in `:uniform` gauge add `chi_i-chi_j` using [`gauge_angles`](@ref).
Theta remains unwrapped. The `Sz_i Sz_j` coupling is never twisted.
"""
function bond_phase(lattice::KagomeCylinder, b::Bond, theta::Real; gauge::Symbol=:seam)
    n = nsites(lattice)
    1 <= b.i <= n && 1 <= b.j <= n ||
        throw(ArgumentError("bond endpoints must belong to the lattice"))
    value = _finite_theta(theta)
    gauge === :seam && return b.wy * value
    gauge === :uniform || throw(ArgumentError("gauge must be :seam or :uniform"))
    si, sj = lattice.sites[b.i], lattice.sites[b.j]
    eta_i = si.y + (si.sublattice === :C ? 0.5 : 0.0)
    eta_j = sj.y + (sj.sublattice === :C ? 0.5 : 0.0)
    return b.wy * value + value * (eta_j - eta_i) / lattice.Ly
end
