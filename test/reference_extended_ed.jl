# Independent geometry/spin-basis reference for the J1-J2-J3 control model.
# Include reference_ed.jl first. No production bond templates or MPOs are used.

function _reference_extended_size(Lx::Integer, Ly::Integer)
    1 <= Lx <= 3 && 3 <= Ly <= 4 ||
        throw(ArgumentError("extended geometry reference supports 1 <= Lx <= 3 and 3 <= Ly <= 4"))
    return Int(Lx), Int(Ly)
end

function _reference_extended_hexagons(Lx::Integer, Ly::Integer)
    Lx, Ly = _reference_extended_size(Lx, Ly)
    # Integer coordinates are (4 Cartesian x, 4 Cartesian y / sqrt(3)).
    # A distance squared is therefore (du^2 + 3dv^2)/16. Two cells of
    # padding contain every NN hexagon touching a retained endpoint, including
    # hexagons across an axial end or a representative periodic seam.
    offsets = ((0, 0), (2, 0), (1, 1))
    vertices = [(; x, y, s, u=4x + 2y + du, v=2y + dv)
                for x in -2:(Lx + 1) for y in -2:(Ly + 1)
                for (s, (du, dv)) in enumerate(offsets)]
    neighbors = [Int[] for _ in vertices]
    for i in eachindex(vertices), j in (i + 1):length(vertices)
        a, b = vertices[i], vertices[j]
        if (a.u - b.u)^2 + 3(a.v - b.v)^2 == 4
            push!(neighbors[i], j)
            push!(neighbors[j], i)
        end
    end
    hexagons = Set{NTuple{6, Int}}()
    function extend_cycle!(path)
        if length(path) == 6
            first(path) in neighbors[last(path)] || return
            # The start is the smallest vertex; fix orientation as well.
            path[2] < path[6] || return
            for a in 1:6, b in (a + 1):6
                (b == a + 1 || (a == 1 && b == 6)) && continue
                path[b] in neighbors[path[a]] && return # a chord
            end
            push!(hexagons, Tuple(path))
            return
        end
        for next in neighbors[last(path)]
            next > first(path) && !(next in path) || continue
            push!(path, next)
            extend_cycle!(path)
            pop!(path)
        end
    end
    for start in eachindex(vertices)
        extend_cycle!([start])
    end
    return (; vertices, hexagons)
end

"""
    reference_extended_bonds(Lx, Ly; J1=1.0, J2=0.5, J3=0.5)

Independent isotropic exchange bonds in canonical `i<j` orientation, including
`wy`, `Jxy`, `Jz`, and `family` (`:J1`, `:J2`, `:J3`). J1/J2 follow Cartesian
distances squared 1/4 and 3/4. J3 joins opposite vertices of chordless NN
hexagons found in an unwrapped plane before endpoint clipping and periodic
identification; distance-one straight-chain neighbors are excluded.

This bounded geometry oracle supports `1≤Lx≤3`, `3≤Ly≤4`. Zero couplings keep
the geometrical bonds, and no finite-cylinder cycle finder or production
templates are used.
"""
function reference_extended_bonds(Lx::Integer, Ly::Integer;
        J1::Real=1.0, J2::Real=0.5, J3::Real=0.5)
    Lx, Ly = _reference_extended_size(Lx, Ly)
    couplings = (; J1=Float64(J1), J2=Float64(J2), J3=Float64(J3))
    all(isfinite, values(couplings)) || throw(ArgumentError("couplings must be finite"))
    keys = Set{Tuple{Int, Int, Int, Symbol}}()
    sites = _reference_sites(Lx, Ly)
    for i in eachindex(sites), j in (i + 1):length(sites), wy in -1:1
        ri, rj = sites[i].position, sites[j].position
        distance2 = (rj[1] + wy * Ly / 2 - ri[1])^2 +
                    (rj[2] + wy * Ly * sqrt(3) / 2 - ri[2])^2
        if isapprox(distance2, 1/4; atol=1e-12, rtol=0)
            push!(keys, (i, j, wy, :J1))
        elseif isapprox(distance2, 3/4; atol=1e-12, rtol=0)
            push!(keys, (i, j, wy, :J2))
        end
    end
    plane = _reference_extended_hexagons(Lx, Ly)
    for cycle in plane.hexagons, k in 1:3
        a, b = plane.vertices[cycle[k]], plane.vertices[cycle[k + 3]]
        0 <= a.x < Lx && 0 <= b.x < Lx || continue
        i = 3(a.x * Ly + mod(a.y, Ly)) + a.s
        j = 3(b.x * Ly + mod(b.y, Ly)) + b.s
        wy = fld(b.y, Ly) - fld(a.y, Ly)
        push!(keys, i < j ? (i, j, wy, :J3) : (j, i, -wy, :J3))
    end
    return [(; i, j, wy, Jxy=getproperty(couplings, family),
             Jz=getproperty(couplings, family), family)
            for (i, j, wy, family) in sort!(collect(keys))]
end

"""
    reference_extended_hamiltonian(Lx, Ly, theta=0.0;
        Q=nothing, nup=nothing, J1=1.0, J2=0.5, J3=0.5,
        gauge=:seam, sparse_matrix=false)

Independent spin-basis matrix with the charge/basis convention of
`reference_hamiltonian`. At most 18 sites and dense dimension 4096 are allowed;
the bounds are checked before geometry or basis construction. Geometry must
also lie within the bounded `reference_extended_bonds` domain. No eigensolve.
"""
function reference_extended_hamiltonian(Lx::Integer, Ly::Integer, theta::Real=0.0;
        Q=nothing, nup=nothing, J1::Real=1.0, J2::Real=0.5, J3::Real=0.5,
        gauge::Symbol=:seam, sparse_matrix::Bool=false)
    N = _reference_site_count(Lx, Ly)
    sector = _reference_sector(N; nup, Q)
    dimension = sector.dimension
    !sparse_matrix && dimension > 4096 &&
        throw(ArgumentError("use sparse_matrix=true for a reference dimension above 4096"))
    isfinite(theta) || throw(ArgumentError("theta must be finite"))
    gauge in (:seam, :uniform) || throw(ArgumentError("gauge must be :seam or :uniform"))
    bonds = reference_extended_bonds(Lx, Ly; J1, J2, J3)
    basis = _reference_basis(N, sector.nup)
    lookup = Dict(state => k for (k, state) in enumerate(basis))
    chi = gauge === :uniform ? _reference_gauge_angles(Lx, Ly, theta) : zeros(N)
    rows, columns, values = Int[], Int[], ComplexF64[]
    # Keep spin-basis arithmetic local so extending the oracle never changes
    # the existing NN reference. A set bit means Up, independently of ITensor.
    for (column, state) in enumerate(basis)
        diagonal = 0.0
        for b in bonds
            zi, zj = _reference_sz_bit(state, b.i), _reference_sz_bit(state, b.j)
            diagonal += b.Jz * zi * zj
            if zi != zj
                flipped = state ⊻ (UInt64(1) << (b.i - 1)) ⊻ (UInt64(1) << (b.j - 1))
                phase = b.wy * theta + chi[b.i] - chi[b.j]
                push!(rows, lookup[flipped])
                push!(columns, column)
                push!(values, (b.Jxy / 2) * cis(zi < zj ? phase : -phase))
            end
        end
        push!(rows, column)
        push!(columns, column)
        push!(values, diagonal)
    end
    H = sparse(rows, columns, values, dimension, dimension)
    return (; H=sparse_matrix ? H : Matrix(H), basis)
end
