# Independent small spin-basis fixtures for an upstream block-sparse
# factorization pathology. All states have total physical Sz = 0; they are
# calibration controls and do not represent the 1/9 kagome target sector.
function _truncation_fixture(probabilities)
    sites = siteinds("S=1/2", 4; conserve_qns=true)
    state = ITensor(ComplexF64, sites)
    configurations = [(1, 2, 1, 2), (1, 1, 2, 2), (2, 2, 1, 1)]
    for (j, (probability, configuration)) in enumerate(zip(probabilities, configurations))
        state[(sites[k] => configuration[k] for k in eachindex(sites))...] =
            sqrt(probability) * cis(0.31j)
    end
    # Three orthogonal left/right configurations have the input Schmidt
    # probabilities exactly. At the central cut their left QN values are
    # 0, +2, -2, so the boundary crosses charge blocks.
    psi = MPS(state, sites; cutoff=0.0, maxdim=4, orthocenter=2)
    return psi
end

@testset "Normalized QN updates preserve reported discarded probability" begin
    cases = ((; name="exact boundary", probabilities=[0.6, 0.2, 0.2], maxdim=2,
                discarded=0.2, decompositions=("svd", "eigen")),
             (; name="near boundary", probabilities=[0.6, 0.2001, 0.1999], maxdim=2,
                discarded=0.1999, decompositions=("svd", "eigen")),
             (; name="resolved boundary", probabilities=[0.6, 0.21, 0.19], maxdim=2,
                discarded=0.19, decompositions=("svd",)),
             (; name="untruncated degeneracy", probabilities=[0.6, 0.2, 0.2], maxdim=3,
                discarded=0.0, decompositions=("eigen",)))
    @testset "$(case.name), $decomposition" for case in cases, decomposition in case.decompositions
        initial = _truncation_fixture(case.probabilities)
        psi = deepcopy(initial)
        phi = psi[2] * psi[3]
        spec = replacebond!(psi, 2, phi; maxdim=case.maxdim, cutoff=0.0,
                            which_decomp=decomposition, ortho="left", normalize=true)
        @test norm(psi) ≈ 1 atol=1e-13
        actual_loss = 1 - abs2(inner(initial, psi))
        observer = KagomeDMRG._SweepDiagnostics()
        # Both kernels must preserve the exact/near-boundary loss after normalization.
        @test length(eigs(spec)) == dim(linkind(psi, 2)) == case.maxdim
        @test actual_loss ≈ case.discarded atol=1e-13
        @test truncerror(spec) ≈ case.discarded atol=1e-13
        @test ITensorMPS.measure!(observer; psi, spec, bond=2, sweep=1,
            half_sweep=1, energy=0.0) === nothing
        @test only(observer.max_truncation_errors) ≈ case.discarded atol=1e-13
    end
end

@testset "Truncation guards still reject invalid backend output" begin
    initial = _truncation_fixture([0.6, 0.2, 0.2])
    _, _, spec = factorize(initial[2] * initial[3], inds(initial[2]);
        maxdim=2, cutoff=0.0, which_decomp="svd")
    # Deliberately pair a rank-two spectrum with a rank-one MPS. This checks
    # the guard independently of whether the installed backend is defective.
    sites = [siteind(initial, j) for j in 1:length(initial)]
    psi = MPS(ComplexF64, sites, ["Up", "Up", "Dn", "Dn"])
    observer = KagomeDMRG._SweepDiagnostics()
    @test_throws ErrorException ITensorMPS.measure!(observer;
        psi, spec, bond=2, sweep=1, half_sweep=1, energy=0.0)
    @test isempty(observer.max_truncation_errors)
    psi[3] *= 0
    @test_throws ErrorException ITensorMPS.measure!(observer;
        psi, spec, bond=2, sweep=1, half_sweep=1, energy=0.0)
    @test isempty(observer.max_truncation_errors)
end

@testset "Truncation rank guard accepts healthy noisy and noiseless DMRG" begin
    # Nonzero noise chooses an eigen-based perturbed density matrix. Its
    # reported error is not asserted to be a wavefunction discarded weight.
    # The noiseless two-spin control is already tested analytically in
    # test_itensor.jl.
    cases = ((4, 0.0), (2, 1e-5), (4, 1e-5))
    @testset "N=$n noise=$noise" for (n, noise) in cases
        (; H, initial) = _heisenberg_chain_control(n; staggered_field=n == 2 ? 3/8 : 0.0)
        observer = KagomeDMRG._SweepDiagnostics()
        _, psi = dmrg(H, initial; nsweeps=2, maxdim=n == 2 ? 1 : 4,
                      cutoff=1e-12, noise, observer, outputlevel=0)
        @test norm(psi) ≈ 1 atol=1e-13
        @test flux(psi) == QN("Sz", 0)
        @test length(observer.max_truncation_errors) == 2
        @test all(isfinite, observer.max_truncation_errors)
        if n == 2
            @test real(inner(psi', H, psi)) ≈ -0.625 atol=1e-13
        else
            @test real(inner(psi', H, psi)) ≈ -(3 + 2sqrt(3)) / 4 atol=1e-11
        end
    end
end
