# Small-system conversions for independent ED checks. These must never be used
# to materialize a research-size MPS or MPO as a full Hilbert-space tensor.
function _contract_small_system(tensors)
    # These bounded ED checks intentionally materialize tensors of high order.
    # Restore the normal warning threshold even if a contraction fails.
    previous = ITensors.disable_warn_order()
    try
        return reduce(*, tensors)
    finally
        ITensors.set_warn_order(previous)
    end
end

function full_mpo_matrix(H, sites)
    length(sites) <= 10 || throw(ArgumentError("dense MPO check is limited to 10 sites"))
    tensor = _contract_small_system(H)
    return reshape(Array(tensor, prime.(sites)..., dag.(sites)...), 2^length(sites), :)
end

function sector_amplitudes(psi, sites, basis)
    length(sites) <= 18 || throw(ArgumentError("dense MPS check is limited to 18 sites"))
    tensor = _contract_small_system(psi)
    amplitudes = vec(Array(tensor, sites...))
    # ITensor Up is index 1; the independent reference encodes Up with bit 1.
    indices = 2^length(sites) .- Int.(basis)
    return amplitudes[indices]
end

function check_small_system(theta; seed=11)
    lattice = kagome_cylinder(1, 3)
    result = run_dmrg(lattice, theta; seed)
    reference = reference_hamiltonian(1, 3, theta)
    eig = eigen(Hermitian(reference.H))
    ground_indices = findall(e -> abs(e - eig.values[1]) < 1e-10, eig.values)
    ground = eig.vectors[:, ground_indices]
    v = sector_amplitudes(result.psi, result.sites, reference.basis)
    projection = ground * (ground' * v)
    leakage = norm(v - projection)
    # At zero flux there are two ground states. Use the ED ground-space
    # projection of the DMRG vector to match the same state for observables.
    matched = projection / norm(projection)
    ed_sz = reference_sz(matched, reference.basis, 9)
    ed_zz = [reference_correlation(matched, reference.basis, 9, i, j)
             for i in 1:9, j in 1:9]
    ed_pm = [reference_correlation(matched, reference.basis, 9, i, j;
                                  operators=("S+", "S-")) for i in 1:9, j in 1:9]
    correlations = spin_correlations(result.psi)
    metrics = (; theta=Float64(theta), seed, ed_energy=eig.values[1],
        dmrg_energy=result.energy, energy_error_per_site=abs(result.energy-eig.values[1])/9,
        ground_degeneracy=length(ground_indices),
        gap_above_ground_space=eig.values[length(ground_indices)+1]-eig.values[1],
        ground_space_leakage=leakage,
        residual_norm=norm(reference.H*v-result.energy*v),
        state_norm_error=abs(norm(v)-1),
        max_sz_error=maximum(abs.(result.sz-ed_sz)),
        max_zz_error=maximum(abs.(correlations.zz-ed_zz)),
        max_pm_error=maximum(abs.(correlations.pm-ed_pm)),
        total_sz=sum(result.sz), variance=result.variance,
        final_sweep_truncation_error=last(result.max_truncation_errors),
        maximum_truncation_error=maximum(result.max_truncation_errors))
    return (; result, reference, metrics)
end
