"""Return local physical Sᶻ expectation values in site order."""
sz_profile(psi::MPS) = real.(expect(psi, "Sz"))

"""
    spin_correlations(psi)

Return matrices `zz[i,j] = <Sz_i Sz_j>` and `pm[i,j] = <S+_i S-_j>`.
The transverse correlations remain complex at nonzero flux.
"""
spin_correlations(psi::MPS) = (; zz=correlation_matrix(psi, "Sz", "Sz"),
                              pm=correlation_matrix(psi, "S+", "S-"))

"""
    bond_energies(psi, lattice, theta=0.0; gauge=:seam)

Return normalized exchange-energy expectations in `lattice.bonds` order.
Each entry is `Jz*<Sz_i Sz_j> + Jxy*real(exp(im*A)*<S+_i S-_j>)` with the
same signed bond phase as the Hamiltonian. Longitudinal fields are excluded:
`sum(bond_energies(...)) - dot(hz, sz_profile(psi))` is the full energy.
Use [`bond_families`](@ref) to group by geometric exchange family.

The input MPS is preserved. This first implementation reuses both complete
correlation matrices, with O(N²) output storage; large-system profiling and
an implementation restricted to the selected bonds remain future work.
"""
function bond_energies(psi::MPS, lattice::KagomeCylinder, theta::Real=0.0;
                       gauge::Symbol=:seam)
    N = nsites(lattice)
    length(psi) == N || throw(ArgumentError("MPS length does not match the lattice"))
    _check_sites([siteind(psi, i) for i in 1:N], N)
    value = _finite_theta(theta)
    gauge in (:seam, :uniform) || throw(ArgumentError("unknown gauge: $gauge"))
    state_norm = norm(psi)
    isfinite(state_norm) && state_norm > 0 ||
        throw(ArgumentError("MPS must have finite positive norm"))
    phases = [bond_phase(lattice, b, value; gauge) for b in lattice.bonds]
    correlations = spin_correlations(psi)
    return [b.Jz * real(correlations.zz[b.i, b.j]) +
            b.Jxy * real(cis(phases[k]) * correlations.pm[b.i, b.j])
            for (k, b) in enumerate(lattice.bonds)]
end

"""
    spin_transfer(lattice, baseline, current; cuts=1:(lattice.Lx-1))

Measure cumulative physical spin transfer relative to a measured zero-flux
profile. At cell cut `c`, the right region has `x >= c`. Return unrounded
left and right changes, and their sum, for every requested cut. This readout
alone does not validate branch continuity or establish a quantized pump.
"""
function spin_transfer(lattice::KagomeCylinder, baseline::AbstractVector{<:Real},
                       current::AbstractVector{<:Real}; cuts=1:(lattice.Lx-1))
    N = nsites(lattice)
    length(baseline) == length(current) == N ||
        throw(ArgumentError("profiles must have one entry per lattice site"))
    all(isfinite, baseline) && all(isfinite, current) ||
        throw(ArgumentError("profiles must be finite"))
    delta = current .- baseline
    return map(collect(cuts)) do cut
        R = right_region(lattice, cut)
        L = setdiff(1:N, R)
        left = sum(delta[L])
        right = sum(delta[R])
        (; cut, left, right, total=left+right)
    end
end
