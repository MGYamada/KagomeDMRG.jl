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

@testset "Finite QN truncation loss is not hidden by normalization" begin
    for decomposition in ("svd", "eigen"),
        probabilities in ([0.6, 0.2, 0.2], [0.6, 0.2001, 0.1999])
        initial = _truncation_fixture(probabilities)
        @test norm(initial) ≈ 1 atol=1e-13
        psi = deepcopy(initial)
        phi = psi[2] * psi[3]
        spec = replacebond!(psi, 2, phi; maxdim=2, cutoff=0.0,
                            which_decomp=decomposition, ortho="left", normalize=true)
        @test norm(psi) ≈ 1 atol=1e-13
        @test norm(psi[3]) > 0
        actual_loss = 1 - abs2(inner(initial, psi))
        observer = KagomeDMRG._SweepDiagnostics()
        if length(eigs(spec)) != dim(linkind(psi, 2))
            # NDTensors 0.4.31 retains only the 0.6 component even though its
            # Spectrum says it retained two components. Normalization and a
            # nonzero-state check alone cannot detect this lost probability.
            @test dim(linkind(psi, 2)) == 1
            @test length(eigs(spec)) == 2
            @test actual_loss ≈ 0.4 atol=1e-13
            @test truncerror(spec) ≈ probabilities[3] atol=1e-13
            @test actual_loss > truncerror(spec) + 0.19
            exception = try
                ITensorMPS.measure!(observer; psi, spec, bond=2, sweep=1,
                    half_sweep=1, energy=0.0)
                nothing
            catch error
                error
            end
            @test exception isa ErrorException
            @test occursin("truncation spectrum", sprint(showerror, exception))
            @test isempty(observer.max_truncation_errors)
        else
            # Permit a future corrected backend when its physical discarded
            # probability agrees with the reported value.
            @test dim(linkind(psi, 2)) == 2
            @test actual_loss ≈ truncerror(spec) atol=1e-13
            @test ITensorMPS.measure!(observer; psi, spec, bond=2, sweep=1,
                half_sweep=1, energy=0.0) === nothing
        end
    end
end

@testset "Resolved QN boundaries preserve reported discarded probability" begin
    for decomposition in ("svd", "eigen"),
        (probabilities, maxdim, discarded) in
            (([0.6, 0.21, 0.19], 2, 0.19), ([0.6, 0.2, 0.2], 3, 0.0))
        initial = _truncation_fixture(probabilities)
        psi = deepcopy(initial)
        spec = replacebond!(psi, 2, psi[2] * psi[3]; maxdim, cutoff=0.0,
                            which_decomp=decomposition, ortho="left", normalize=true)
        @test length(eigs(spec)) == dim(linkind(psi, 2)) == maxdim
        @test 1 - abs2(inner(initial, psi)) ≈ discarded atol=1e-13
        @test truncerror(spec) ≈ discarded atol=1e-13
        observer = KagomeDMRG._SweepDiagnostics()
        @test ITensorMPS.measure!(observer; psi, spec, bond=2, sweep=1,
            half_sweep=1, energy=0.0) === nothing
        @test only(observer.max_truncation_errors) ≈ discarded atol=1e-13
    end
end

@testset "Truncation rank guard accepts healthy noisy and noiseless DMRG" begin
    # Nonzero noise chooses an eigen-based perturbed density matrix. Its
    # reported error is not asserted to be a wavefunction discarded weight.
    for noise in (0.0, 1e-5), n in (2, 4)
        sites = siteinds("S=1/2", n; conserve_qns=true)
        terms = OpSum()
        for j in 1:(n - 1)
            terms += "Sz", j, "Sz", j + 1
            terms += 0.5, "S+", j, "S-", j + 1
            terms += 0.5, "S-", j, "S+", j + 1
        end
        if n == 2
            terms += -3 / 8, "Sz", 1
            terms += 3 / 8, "Sz", 2
        end
        H = MPO(ComplexF64, terms, sites)
        initial = MPS(ComplexF64, sites, [isodd(j) ? "Up" : "Dn" for j in 1:n])
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
