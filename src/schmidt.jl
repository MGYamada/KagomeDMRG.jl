"""
    schmidt_diagnostics(psi::MPS, bond::Integer)

Measure the Schmidt probabilities and absolute left-block integer charges
`left_q = 2Sz_left` at the MPS bond after site `bond` (`1 <= bond < N`).
The left region is the sites `1:bond`; a geometric cylinder cut `c` corresponds
to `bond = 3lattice.Ly*c` for the documented site ordering.

Return `bond`, `probabilities` (descending, normalized across all sectors),
`left_q` (one charge per probability), the natural-log `entropy`,
`mean_left_sz`, `variance_left_sz`, and the input `norm_squared`.
Zero singular values contribute zero to the entropy. No sector multiplicity
is applied: every singular value already represents one Schmidt state.

The input must be an open-boundary QN spin-1/2 MPS with outgoing physical
indices in the `Sz = 2Sᶻ` convention. A private copy is canonicalized and
decomposed without a cutoff or dimension limit; the input is unchanged and
no full Hilbert-space tensor is constructed. Link QN signs and additive
offsets are removed using tensor flux conservation, not a fit to measured
spin densities. This diagnoses a state; it does not certify continuation or
quantized transport.
"""
function schmidt_diagnostics(psi::MPS, bond::Integer)
    N = length(psi)
    1 <= bond < N || throw(ArgumentError("bond must lie in 1:(length(psi)-1)"))
    b = Int(bond)
    sites = [siteind(psi, j) for j in 1:N]
    all(s -> s isa Index, sites) ||
        throw(ArgumentError("expected one physical index per MPS tensor"))
    _check_sites(sites, N)
    all(s -> dir(s) == ITensors.Out, sites) ||
        throw(ArgumentError("physical spin indices must point out of the MPS"))
    for j in 1:N
        order(psi[j]) == (j in (1, N) ? 2 : 3) ||
            throw(ArgumentError("expected an open MPS with one link per neighboring pair"))
        tensor_norm = norm(psi[j])
        isfinite(tensor_norm) && tensor_norm > 0 ||
            throw(ArgumentError("MPS tensors must have finite nonzero norm"))
        ITensors.checkflux(psi[j])
    end
    all(j -> length(commoninds(psi[j], psi[j+1])) == 1, 1:(N-1)) ||
        throw(ArgumentError("expected one link per neighboring MPS pair"))

    state = deepcopy(psi)
    orthogonalize!(state, b)
    center_norm = norm(state[b])
    isfinite(center_norm) && center_norm > 0 ||
        throw(ArgumentError("MPS must have finite nonzero norm"))

    # Omitting both cutoff and maxdim selects the upstream no-truncation
    # branch, including at degenerate QN-sector boundaries.
    left_indices = uniqueinds(state[b], state[b+1])
    U, S, V = svd(state[b], left_indices)
    left_link = commonind(U, S)
    right_link = commonind(V, S)
    weights = [abs2(S[left_link => a, right_link => a]) for a in 1:dim(left_link)]
    norm_squared = sum(weights)
    isfinite(norm_squared) && norm_squared > 0 ||
        throw(ArgumentError("Schmidt decomposition has finite nonzero norm requirement"))

    # Sum local tensor fluxes over the left isometry. Interior link charges
    # cancel in pairs. With outgoing physical sites, conservation gives
    # physical q_left + directed q_boundary = total tensor flux on the left.
    # The offset is needed even when the total MPS charge is nonzero or its
    # charge has been redistributed among tensors by a link relabeling.
    left_flux = flux(U)
    for j in 1:(b-1)
        left_flux += flux(state[j])
    end
    left_q = [val(left_flux - flux(left_link => a), "Sz") for a in 1:dim(left_link)]
    all(q -> abs(q) <= b && iseven(q-b), left_q) ||
        error("Schmidt charges are inconsistent with a block of $b spin-1/2 sites")

    permutation = sortperm(weights; rev=true)
    probabilities = weights[permutation] ./ norm_squared
    left_q = left_q[permutation]
    mean_left_sz = sum(probabilities .* left_q) / 2
    variance_left_sz = sum(probabilities .* (left_q ./ 2 .- mean_left_sz).^2)
    entropy = -sum(p -> iszero(p) ? zero(p) : p*log(p), probabilities)
    return (; bond=b, probabilities, left_q, entropy, mean_left_sz,
            variance_left_sz, norm_squared)
end
