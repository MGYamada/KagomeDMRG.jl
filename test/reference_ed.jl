# Independent spin-basis reference for small-system tests. This file deliberately
# does not import KagomeDMRG or use its lattice templates, phases, or MPO builder.
using LinearAlgebra
using SparseArrays

function _reference_sites(Lx::Integer, Ly::Integer)
    Lx >= 1 || throw(ArgumentError("Lx must be positive"))
    Ly >= 3 || throw(ArgumentError("the reference requires Ly >= 3"))
    sites = NamedTuple[]
    offsets = ((0.0, 0.0), (0.5, 0.0), (0.25, sqrt(3) / 4))
    for x in 0:(Lx - 1), y in 0:(Ly - 1), (s, offset) in enumerate(offsets)
        push!(sites, (; x, y, sublattice=(:A, :B, :C)[s],
            position=(x + y / 2 + offset[1], sqrt(3) * y / 2 + offset[2])))
    end
    return sites
end

"""
    reference_bonds(Lx, Ly; Jxy=1.0, Jz=1.0)

Find nearest neighbors by Cartesian distance over periodic images `-1:1`.
The canonical orientation is `i < j`; `wy` shifts the image of site `j` by
`wy * Ly * a2`. No production bond templates are used.
"""
function reference_bonds(Lx::Integer, Ly::Integer; Jxy::Real=1.0, Jz::Real=1.0)
    sites = _reference_sites(Lx, Ly)
    bonds = NamedTuple{(:i, :j, :wy, :Jxy, :Jz), Tuple{Int, Int, Int, Float64, Float64}}[]
    for i in eachindex(sites), j in (i + 1):length(sites), wy in -1:1
        ri, rj = sites[i].position, sites[j].position
        dx = rj[1] + wy * Ly / 2 - ri[1]
        dy = rj[2] + wy * Ly * sqrt(3) / 2 - ri[2]
        if isapprox(dx^2 + dy^2, 1 / 4; atol=1e-12, rtol=0)
            push!(bonds, (; i, j, wy, Jxy=Float64(Jxy), Jz=Float64(Jz)))
        end
    end
    return bonds
end

function _reference_basis(N::Integer, nup)
    1 <= N <= 18 || throw(ArgumentError("this test reference supports 1 <= N <= 18"))
    if isnothing(nup)
        N % 9 == 0 || throw(ArgumentError("the 1/9 target requires N divisible by 9"))
        nup = 5 * (N ÷ 9)
    end
    if nup === :all
        return collect(UInt64(0):((UInt64(1) << N) - 1))
    end
    nup isa Integer || throw(ArgumentError("nup must be an integer, :all, or nothing"))
    0 <= nup <= N || throw(ArgumentError("nup must lie between zero and N"))
    return [state for state in UInt64(0):((UInt64(1) << N) - 1) if count_ones(state) == nup]
end

_reference_sz_bit(state::UInt64, i::Integer) = iszero(state & (UInt64(1) << (i - 1))) ? -0.5 : 0.5

function _reference_gauge_angles(Lx::Integer, Ly::Integer, theta::Real)
    return [-theta * (s.y + (s.sublattice === :C ? 0.5 : 0.0)) / Ly
            for s in _reference_sites(Lx, Ly)]
end

"""
    reference_hamiltonian(Lx, Ly, theta=0.0;
                          nup=nothing, Jxy=1.0, Jz=1.0,
                          gauge=:seam, sparse_matrix=false)

Return `(H, basis)` as a named tuple. Bit `i-1` is one for an Up spin at site
`i`, and basis integers are increasing. `nup=nothing` selects the 1/9 target;
`nup=:all` constructs all magnetization sectors for explicit U(1) checks.
This reference supports at most 18 sites and refuses dense matrices above
dimension 4096. It performs no eigensolve.
"""
function reference_hamiltonian(Lx::Integer, Ly::Integer, theta::Real=0.0;
        nup=nothing, Jxy::Real=1.0, Jz::Real=1.0, gauge::Symbol=:seam,
        sparse_matrix::Bool=false)
    gauge in (:seam, :uniform) || throw(ArgumentError("gauge must be :seam or :uniform"))
    N = 3 * Lx * Ly
    bonds = reference_bonds(Lx, Ly; Jxy, Jz)
    basis = _reference_basis(N, nup)
    dimension = length(basis)
    !sparse_matrix && dimension > 4096 &&
        throw(ArgumentError("use sparse_matrix=true for a reference dimension above 4096"))
    lookup = Dict(state => k for (k, state) in enumerate(basis))
    chi = gauge === :uniform ? _reference_gauge_angles(Lx, Ly, theta) : zeros(N)
    rows, columns, values = Int[], Int[], ComplexF64[]
    for (column, state) in enumerate(basis)
        diagonal = 0.0
        for b in bonds
            zi, zj = _reference_sz_bit(state, b.i), _reference_sz_bit(state, b.j)
            diagonal += b.Jz * zi * zj
            if zi != zj
                flipped = state ⊻ (UInt64(1) << (b.i - 1)) ⊻ (UInt64(1) << (b.j - 1))
                # S+_i S-_j acts on (down,up); its adjoint acts on (up,down).
                phase = b.wy * theta + chi[b.i] - chi[b.j]
                amplitude = (b.Jxy / 2) * cis(zi < zj ? phase : -phase)
                push!(rows, lookup[flipped])
                push!(columns, column)
                push!(values, amplitude)
            end
        end
        push!(rows, column)
        push!(columns, column)
        push!(values, diagonal)
    end
    H = sparse(rows, columns, values, dimension, dimension)
    return (; H=sparse_matrix ? H : Matrix(H), basis)
end

"""Diagonal entries of `U = exp(i sum(chi_i Sz_i))` in the given basis."""
function reference_gauge_diagonal(basis, Lx::Integer, Ly::Integer, theta::Real)
    chi = _reference_gauge_angles(Lx, Ly, theta)
    return [cis(sum(chi[i] * _reference_sz_bit(state, i) for i in eachindex(chi)))
            for state in basis]
end

function _reference_check_state(psi, basis, N::Integer)
    length(psi) == length(basis) || throw(DimensionMismatch("state and basis lengths differ"))
    1 <= N <= 18 || throw(ArgumentError("this test reference supports 1 <= N <= 18"))
    all(s -> s < UInt64(1) << N, basis) || throw(ArgumentError("basis state exceeds N sites"))
    norm_squared = sum(abs2, psi)
    norm_squared > 0 || throw(ArgumentError("state must have nonzero norm"))
    return norm_squared
end

"""Normalized physical `Sz` expectation at every site."""
function reference_sz(psi, basis, N::Integer)
    norm_squared = _reference_check_state(psi, basis, N)
    return [sum(abs2(psi[k]) * _reference_sz_bit(state, i)
                for (k, state) in enumerate(basis)) / norm_squared for i in 1:N]
end

function _reference_apply_spin(state::UInt64, site::Integer, op)
    op = String(op)
    z = _reference_sz_bit(state, site)
    op == "Sz" && return (state, z)
    op in ("Id", "I") && return (state, 1.0)
    if op == "S+"
        return z < 0 ? (state ⊻ (UInt64(1) << (site - 1)), 1.0) : (state, 0.0)
    elseif op == "S-"
        return z > 0 ? (state ⊻ (UInt64(1) << (site - 1)), 1.0) : (state, 0.0)
    end
    throw(ArgumentError("unsupported reference operator $op"))
end

"""
    reference_correlation(psi, basis, N, i, j; operators=("Sz", "Sz"))

Normalized expectation of `operators[1]_i * operators[2]_j`. Operators act
right to left even if `i == j`; supported names are `Sz`, `S+`, `S-`, `Id`.
"""
function reference_correlation(psi, basis, N::Integer, i::Integer, j::Integer;
        operators=("Sz", "Sz"))
    norm_squared = _reference_check_state(psi, basis, N)
    1 <= i <= N && 1 <= j <= N || throw(ArgumentError("site index lies outside the system"))
    length(operators) == 2 || throw(ArgumentError("two operator names are required"))
    lookup = Dict(state => k for (k, state) in enumerate(basis))
    expectation = 0.0 + 0.0im
    for (column, state) in enumerate(basis)
        intermediate, a = _reference_apply_spin(state, j, operators[2])
        iszero(a) && continue
        final, b = _reference_apply_spin(intermediate, i, operators[1])
        iszero(b) && continue
        row = get(lookup, final, 0)
        iszero(row) && continue
        expectation += conj(psi[row]) * a * b * psi[column]
    end
    return expectation / norm_squared
end
