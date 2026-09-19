"""
    KagomeTriangle(sites, image_y, kind)

Three ordered, distinct positive site indices and their periodic images in
units of the wrap vector `Ly*a2`. `kind` is `:up` or `:down`. The constructor
validates the data, not geometric orientation; [`oriented_triangles`](@ref)
supplies counterclockwise (CCW) elementary triangles on the canonical cylinder.
"""
struct KagomeTriangle
    sites::NTuple{3,Int}
    image_y::NTuple{3,Int}
    kind::Symbol

    function KagomeTriangle(sites::NTuple{3,<:Integer}, image_y::NTuple{3,<:Integer},
                            kind::Symbol)
        all(i -> !(i isa Bool) && 0 < i <= typemax(Int), sites) &&
            length(unique(sites)) == 3 ||
            throw(ArgumentError("triangle sites must be distinct positive integers fitting in Int"))
        all(w -> !(w isa Bool) && typemin(Int) <= w <= typemax(Int), image_y) ||
            throw(ArgumentError("triangle images must be integers fitting in Int"))
        kind in (:up, :down) || throw(ArgumentError("triangle kind must be :up or :down"))
        return new(Int.(sites), Int.(image_y), kind)
    end
end

"""
    oriented_triangles(lattice) -> Vector{KagomeTriangle}

Return the `(2Lx-1)Ly` elementary triangles with CCW orientation in the
unwrapped plane. In increasing `x`, then `y`, emit the up triangle
`A(x,y), B(x,y), C(x,y)`, followed when `x>=1` by the down triangle
`A(x,y), B(x-1,y), C(x,y-1)`. The latter C vertex has image `-1` at `y=0`.

Require the canonical sites and ordering of [`kagome_cylinder`](@ref).
The list is geometric: changing exchange couplings or adding J2/J3 bonds
does not add triangles or alter their orientation.
"""
function oriented_triangles(lattice::KagomeCylinder)::Vector{KagomeTriangle}
    reference = kagome_cylinder(lattice.Lx, lattice.Ly)
    lattice.sites == reference.sites ||
        throw(ArgumentError("triangles require the documented kagome sites and ordering"))
    triangles = KagomeTriangle[]
    sizehint!(triangles, (2lattice.Lx - 1) * lattice.Ly)
    for x in 0:(lattice.Lx - 1), y in 0:(lattice.Ly - 1)
        a = site_index(lattice, x, y, :A)
        push!(triangles, KagomeTriangle((a, site_index(lattice, x, y, :B),
            site_index(lattice, x, y, :C)), (0, 0, 0), :up))
        if x >= 1
            image, wrapped_y = fldmod(y - 1, lattice.Ly)
            push!(triangles, KagomeTriangle((a, site_index(lattice, x - 1, y, :B),
                site_index(lattice, x, wrapped_y, :C)), (0, 0, image), :down))
        end
    end
    return triangles
end

"""
    scalar_chirality_mpo(sites, vertices; angles=(0.0, 0.0, 0.0))

Return a `ComplexF64` MPO for `U * Si⋅(Sj×Sk) * U†`, with `(i,j,k)=vertices`
and `U=exp(im*sum(angles[a]*Sz_vertices[a]))`. Exchanging two vertices reverses
the undressed operator's sign. The six ladder terms conserve total integer
charge `Q=2Sz`; each cyclic term is `(im/2) Sz_a (S+_b S-_c - S-_b S+_c)`.
Angles follow the supplied vertex order, and rotate the plus/minus term by
`cis(angles[b]-angles[c])`. No term is added to the model Hamiltonian.
"""
function scalar_chirality_mpo(sites, vertices::NTuple{3,<:Integer};
                              angles=(0.0, 0.0, 0.0))
    N = length(sites)
    _check_sites(sites, N)
    all(i -> !(i isa Bool) && 1 <= i <= N, vertices) &&
        length(unique(vertices)) == 3 ||
        throw(ArgumentError("chirality requires three distinct valid site indices"))
    (angles isa Tuple || angles isa AbstractVector) && length(angles) == 3 &&
        all(a -> a isa Real, angles) ||
        throw(ArgumentError("chirality angles must contain three finite real values"))
    alpha = Float64.(angles)
    all(isfinite, alpha) || throw(ArgumentError("chirality angles must be finite"))
    terms = OpSum()
    for (a, b, c) in ((1, 2, 3), (2, 3, 1), (3, 1, 2))
        difference = alpha[b] - alpha[c]
        isfinite(difference) || throw(ArgumentError("chirality angle difference must be finite"))
        coefficient = 0.5im * cis(difference)
        terms += coefficient, "Sz", vertices[a], "S+", vertices[b], "S-", vertices[c]
        terms += conj(coefficient), "Sz", vertices[a], "S-", vertices[b], "S+", vertices[c]
    end
    return MPO(ComplexF64, terms, sites)
end

"""
    triangle_chiralities(psi, lattice, theta=0.0; gauge=:seam)

Return normalized real chirality expectations in [`oriented_triangles`](@ref)
order, preserving the input MPS. At zero flux these are ordinary CCW scalar
chiralities. At finite flux they are locally transported (dressed) operators:
`alpha[a]=-image_y[a]*theta`, plus `gauge_angles(lattice,theta)[sites[a]]`
in uniform gauge. Thus a plus/minus pair has the same signed Peierls phase
as the corresponding oriented exchange bond. Gauge comparisons require the
state to be rotated with the same `U`; the bare finite-flux operator is different.

This reference readout builds and contracts one small MPO per triangle on a
private normalized MPS copy. It neither constructs a full Hilbert-space tensor
nor certifies a chiral phase; large-system measurement cost is not calibrated.
"""
function triangle_chiralities(psi::MPS, lattice::KagomeCylinder, theta::Real=0.0;
                               gauge::Symbol=:seam)
    N = nsites(lattice)
    length(psi) == N || throw(ArgumentError("MPS length does not match the lattice"))
    sites = [siteind(psi, i) for i in 1:N]
    _check_sites(sites, N)
    value = _finite_theta(theta)
    gauge in (:seam, :uniform) || throw(ArgumentError("unknown gauge: $gauge"))
    triangles = oriented_triangles(lattice)
    rotations = gauge === :uniform ? gauge_angles(lattice, value) : zeros(N)
    state = deepcopy(psi)
    state_norm = norm(state)
    isfinite(state_norm) && state_norm > 0 ||
        throw(ArgumentError("MPS must have finite positive norm"))
    normalize!(state)
    normalization = real(inner(state, state))
    values = Float64[]
    sizehint!(values, length(triangles))
    for triangle in triangles
        angles = ntuple(a -> -triangle.image_y[a] * value + rotations[triangle.sites[a]], 3)
        operator = scalar_chirality_mpo(sites, triangle.sites; angles)
        measured = inner(state', operator, state) / normalization
        isfinite(measured) && abs(imag(measured)) <= 1e-10 * max(1, abs(real(measured))) ||
            error("chirality expectation is nonfinite or has a significant imaginary part")
        push!(values, real(measured))
    end
    return values
end
