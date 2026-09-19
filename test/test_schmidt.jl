# These dense helpers are independent spin-basis checks, restricted to N <= 9.
# Bits encode Up = 1, and bit j-1 is physical site j.
function _schmidt_test_mps(amplitudes, basis, sites)
    N = length(sites)
    N <= 9 || throw(ArgumentError("Schmidt reference is limited to 9 sites"))
    tensor = ITensor(ComplexF64, sites...)
    for (amplitude, bits) in zip(amplitudes, basis)
        indices = [sites[j] => (isodd(bits >> (j-1)) ? 1 : 2) for j in 1:N]
        tensor[indices...] = amplitude
    end
    return MPS(tensor, sites)
end

function _spin_basis_schmidt(amplitudes, basis, N, bond)
    N <= 9 || throw(ArgumentError("Schmidt reference is limited to 9 sites"))
    matrix = zeros(ComplexF64, 2^bond, 2^(N-bond))
    for (amplitude, bits) in zip(amplitudes, basis)
        left = Int(bits & (2^bond-1)) + 1
        right = Int(bits >> bond) + 1
        matrix[left, right] = amplitude
    end
    matrix ./= norm(matrix)
    density = matrix * matrix'
    charges = [2count_ones(left)-bond for left in 0:(2^bond-1)]
    probabilities, left_q = Float64[], Int[]
    for q in sort(unique(charges))
        rows = findall(==(q), charges)
        weights = abs2.(svdvals(matrix[rows, :]))
        append!(probabilities, weights)
        append!(left_q, fill(q, length(weights)))
    end
    return (; probabilities, left_q,
        mean_left_sz=real(sum(diag(density) .* charges))/2,
        entropy=-sum(p -> iszero(p) ? 0.0 : p*log(p), probabilities))
end

function _relabel_test_link(psi, bond; shift, reverse_arrow)
    # Change the QN labels and arrow without changing any tensor entries.
    # Dense storage is used only to construct this bounded adversarial fixture.
    length(psi) <= 9 || throw(ArgumentError("link relabel fixture is limited to 9 sites"))
    result = deepcopy(psi)
    old = commonind(result[bond], result[bond+1])
    sign = reverse_arrow ? -1 : 1
    spaces = [QN("Sz", sign*val(q, "Sz")+shift) => d for (q, d) in space(old)]
    replacement = Index(spaces; tags=tags(old),
                        dir=reverse_arrow ? dir(dag(old)) : dir(old))
    for (site, link) in ((bond, replacement), (bond+1, dag(replacement)))
        original_indices = collect(inds(result[site]))
        new_indices = [index == old ? link : index for index in original_indices]
        result[site] = ITensor(Array(result[site], original_indices...), new_indices...)
    end
    return result
end

@testset "Absolute Schmidt charge of product states" begin
    sites = siteinds("S=1/2", 9; conserve_qns=true)
    cases = (("mixed", [fill("Up", 5); fill("Dn", 4)], 1, 4),
             ("alternating", ["Dn", "Up", "Dn", "Up", "Dn", "Up", "Dn", "Up", "Up"], 9, 8),
             ("all up", fill("Up", 9), 5, 1), ("all down", fill("Dn", 9), 1, 8))
    @testset "$name" for (name, labels, center, bond) in cases
        psi = MPS(ComplexF64, sites, labels)
        orthogonalize!(psi, center)
        before = deepcopy(psi)
        limits = (ITensorMPS.leftlim(psi), ITensorMPS.rightlim(psi))
        result = schmidt_diagnostics(psi, bond)
        q = sum(label == "Up" ? 1 : -1 for label in labels[1:bond])
        @test result.probabilities ≈ [1.0] atol=1e-14
        @test result.left_q == [q]
        @test result.mean_left_sz ≈ q/2 atol=1e-14
        @test result.variance_left_sz ≈ 0 atol=1e-14
        @test result.entropy ≈ 0 atol=1e-14
        @test (ITensorMPS.leftlim(psi), ITensorMPS.rightlim(psi)) == limits
        @test all(inds(psi[j]) == inds(before[j]) &&
                  norm(psi[j]-before[j]) == 0 for j in 1:9)
    end
end

@testset "Known complex Schmidt weights and input validation" begin
    sites = siteinds("S=1/2", 2; conserve_qns=true)
    for p in (0.8, 0.5, 1.0)
        psi = _schmidt_test_mps([sqrt(p), cis(0.63)*sqrt(1-p)], [1, 2], sites)
        before = deepcopy(psi)
        result = schmidt_diagnostics(psi, 1)
        @test sum(result.probabilities) ≈ 1 atol=1e-14
        @test sum(result.probabilities[result.left_q .== 1]) ≈ p atol=1e-14
        @test sum(result.probabilities[result.left_q .== -1]) ≈ 1-p atol=1e-14
        @test result.mean_left_sz ≈ p-0.5 atol=1e-14
        @test result.variance_left_sz ≈ p*(1-p) atol=1e-14
        expected_entropy = p == 1 ? 0.0 : -p*log(p)-(1-p)*log(1-p)
        @test result.entropy ≈ expected_entropy atol=1e-14
        @test all(norm(psi[j]-before[j]) == 0 for j in 1:2)
        if p == 0.8
            scaled = deepcopy(psi)
            scaled[1] *= 2cis(0.37)
            scaled_result = schmidt_diagnostics(scaled, 1)
            @test scaled_result.norm_squared ≈ 4 atol=1e-13
            @test scaled_result.probabilities ≈ result.probabilities atol=1e-14
            @test scaled_result.left_q == result.left_q
        end
    end
    product = MPS(ComplexF64, sites, ["Up", "Dn"])
    @test_throws ArgumentError schmidt_diagnostics(product, 0)
    @test_throws ArgumentError schmidt_diagnostics(product, 2)
    @test_throws ArgumentError schmidt_diagnostics(MPS(ComplexF64,
        siteinds("S=1/2", 2), ["Up", "Dn"]), 1)
    zero_state = deepcopy(product)
    zero_state[1] *= 0
    @test_throws ArgumentError schmidt_diagnostics(zero_state, 1)
    nonfinite = deepcopy(product)
    nonfinite[1] *= NaN
    @test_throws ArgumentError schmidt_diagnostics(nonfinite, 1)
end

@testset "9-site charge-resolved Schmidt spectrum against spin-basis SVD" begin
    sites = siteinds("S=1/2", 9; conserve_qns=true)
    let theta = 0.37
        reference = reference_hamiltonian(1, 3, theta)
        eigenstates = eigen(Hermitian(reference.H)).vectors
        # A complex superposition tests the full fixed-Q sector, independently
        # of DMRG convergence and of the degeneracy of the zero-flux ground state.
        amplitudes = (eigenstates[:, 1] + (0.31+0.27im)*eigenstates[:, 3])
        amplitudes ./= norm(amplitudes)
        psi = _schmidt_test_mps(amplitudes, reference.basis, sites)
        @test flux(psi) == QN("Sz", 1)
        densities = reference_sz(amplitudes, reference.basis, 9)
        # Center left of the cut, on the cut, and right of the cut.
        @testset "center=$center bond=$bond" for (center, bond) in ((1, 8), (4, 4), (9, 1))
            orthogonalize!(psi, center)
            result = schmidt_diagnostics(psi, bond)
            exact = _spin_basis_schmidt(amplitudes, reference.basis, 9, bond)
            @test result.norm_squared ≈ 1 atol=1e-12
            @test sum(result.probabilities) ≈ 1 atol=1e-13
            @test issorted(result.probabilities; rev=true)
            @test result.mean_left_sz ≈ exact.mean_left_sz atol=1e-12
            @test result.mean_left_sz ≈ sum(densities[1:bond]) atol=1e-12
            @test result.entropy ≈ exact.entropy atol=1e-12
            for q in unique(exact.left_q)
                actual = sort(filter(>(1e-14), result.probabilities[result.left_q .== q]))
                expected = sort(filter(>(1e-14), exact.probabilities[exact.left_q .== q]))
                @test actual ≈ expected atol=1e-12
            end
        end
        # Explicitly shift a link's origin and reverse its arrow. Absolute
        # left charge must remain calibrated without using the density profile.
        bonds = (3, 4, 5)
        baseline = Dict(b => schmidt_diagnostics(psi, b) for b in bonds)
        for reverse_arrow in (false, true)
            relabeled = _relabel_test_link(psi, 4; shift=7, reverse_arrow)
            @test abs(inner(relabeled, psi)) ≈ 1 atol=1e-12
            for bond in bonds
                result = schmidt_diagnostics(relabeled, bond)
                @test result.probabilities ≈ baseline[bond].probabilities atol=1e-12
                @test result.left_q == baseline[bond].left_q
                @test result.mean_left_sz ≈ baseline[bond].mean_left_sz atol=1e-12
            end
        end
    end
end

@testset "Geometric spin transfer agrees with Schmidt charge on two cuts" begin
    # Synthetic readout calibration, not a flux-evolved state or a quantized
    # pump claim. Both configurations have Nup=15, Ndown=12 and Q=3.
    # One Up spin moves from the left end to the right end across both cuts.
    lattice = kagome_cylinder(3, 3)
    sites = spin_sites(lattice)
    labels = [fill("Up", 15); fill("Dn", 12)]
    moved_labels = copy(labels)
    moved_labels[1], moved_labels[end] = "Dn", "Up"
    baseline = MPS(ComplexF64, sites, labels)
    moved = MPS(ComplexF64, sites, moved_labels)
    baseline_profile = sz_profile(baseline)
    @test flux(baseline) == flux(moved) == QN("Sz", 3)
    @test sum(baseline_profile) ≈ 1.5 atol=1e-13
    baseline_schmidt = [schmidt_diagnostics(baseline, 3lattice.Ly*cut) for cut in 1:2]

    for p in (0.2, 1.0)
        # Direct-sum addition preserves both branches exactly at bond dimension
        # two. It does not materialize the 27-site Hilbert space or run DMRG.
        current = if p == 1
            deepcopy(moved)
        else
            add(sqrt(1-p)*baseline, cis(0.63)*sqrt(p)*moved; alg="directsum")
        end
        @test maxlinkdim(current) <= 2
        @test flux(current) == QN("Sz", 3)
        @test norm(current) ≈ 1 atol=1e-13
        current_profile = sz_profile(current)
        expected_delta = zeros(27)
        expected_delta[1], expected_delta[end] = -p, p
        @test current_profile-baseline_profile ≈ expected_delta atol=1e-13
        @test sum(current_profile) ≈ sum(baseline_profile) atol=1e-13
        transfer = spin_transfer(lattice, baseline_profile, current_profile)
        @test length(transfer) == 2
        expected_entropy = p in (0, 1) ? 0.0 : -p*log(p)-(1-p)*log(1-p)
        for cut in 1:2
            bond = 3lattice.Ly*cut
            @test right_region(lattice, cut) == collect((bond+1):27)
            diag = schmidt_diagnostics(current, bond)
            delta_left = diag.mean_left_sz-baseline_schmidt[cut].mean_left_sz
            @test diag.mean_left_sz ≈ sum(current_profile[1:bond]) atol=1e-13
            @test diag.variance_left_sz ≈ p*(1-p) atol=1e-13
            @test diag.entropy ≈ expected_entropy atol=1e-13
            @test delta_left ≈ -p atol=1e-13
            @test transfer[cut].left ≈ delta_left atol=1e-13
            @test transfer[cut].right ≈ p atol=1e-13
            @test transfer[cut].total ≈ 0 atol=1e-13
        end
    end
end
