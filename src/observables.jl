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
