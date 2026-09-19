@testset "Complex U(1) MPO versus independent ED" begin
    lattice = kagome_cylinder(1, 3)
    sites = spin_sites(lattice)
    @test all(hasqns, sites)
    cases = ((0.0, :seam), (0.37, :seam), (0.37, :uniform), (-0.71, :uniform))
    seam = nothing
    @testset "theta=$theta gauge=$gauge" for (theta, gauge) in cases
        H = twisted_exchange_mpo(sites, lattice, theta; gauge)
        @test all(eltype(t) == ComplexF64 for t in H)
        full = full_mpo_matrix(H, sites)
        theta == 0.37 && gauge === :seam && (seam = full)
        reference = reference_hamiltonian(1, 3, theta; gauge, nup=:all)
        indices = 2^9 .- Int.(reference.basis)
        reordered = full[indices, indices]
        @test reordered ≈ reference.H atol=2e-12 rtol=0
        @test full ≈ full' atol=2e-12 rtol=0
        Q = Diagonal([2count_ones(s)-9 for s in reference.basis])
        @test norm(Q*reordered-reordered*Q) < 1e-12
    end
    periodic = full_mpo_matrix(twisted_exchange_mpo(sites, lattice, 0.37+2pi), sites)
    @test seam ≈ periodic atol=2e-12 rtol=0

    # Generic XXZ coefficients expose missing factors even away from SU(2).
    xxz = kagome_cylinder(1, 3; Jxy=0.7, Jz=1.3)
    Hxxz = full_mpo_matrix(twisted_exchange_mpo(sites, xxz, 0.37), sites)
    reference = reference_hamiltonian(1, 3, 0.37; Jxy=0.7, Jz=1.3, nup=:all)
    indices = 2^9 .- Int.(reference.basis)
    @test Hxxz[indices, indices] ≈ reference.H atol=2e-12 rtol=0

    @test_throws ArgumentError twisted_exchange_mpo(sites, lattice, NaN)
    @test_throws ArgumentError twisted_exchange_mpo(sites, lattice, 0; gauge=:invalid)
    @test_throws ArgumentError twisted_exchange_mpo(sites, lattice, 0; hz=zeros(8))
    @test_throws ArgumentError twisted_exchange_mpo(siteinds("S=1/2", 9), lattice, 0)
end

@testset "Fixed-charge complex initialization and warm-start contract" begin
    lattice = kagome_cylinder(1, 3)
    sites = spin_sites(lattice)
    p = initial_mps(sites; seed=11)
    same = initial_mps(sites; seed=11)
    other = initial_mps(sites; seed=29)
    @test flux(p) == QN("Sz", 1)
    @test all(eltype(t) == ComplexF64 for t in p)
    @test abs(inner(p, same)) ≈ 1 atol=1e-12
    @test abs(inner(p, other)) < 1-1e-5
    @test sum(sz_profile(p)) ≈ 0.5 atol=1e-12
    zero = initial_mps(spin_sites(kagome_cylinder(1, 4)); Q=0, seed=11)
    @test flux(zero) == QN("Sz", 0)
    @test sum(sz_profile(zero)) ≈ 0 atol=1e-12
    @test_throws ArgumentError initial_mps(sites; linkdim=0)
    @test_throws ArgumentError run_dmrg(lattice, 0; sites, nsweeps=0)
    @test_throws ArgumentError run_dmrg(lattice, 0; sites, maxdim=Int[])
    @test_throws ArgumentError run_dmrg(lattice, 0; sites, cutoff=-1.0)
    @test_throws ArgumentError run_dmrg(lattice, 0; sites, psi0=initial_mps(spin_sites(lattice)))
    wrong_charge = MPS(ComplexF64, sites, fill("Up", 9))
    @test_throws ArgumentError run_dmrg(lattice, 0; sites, psi0=wrong_charge)
    @test_throws ArgumentError run_dmrg(lattice, 0; sites, psi0=p, Q=-1)
end

@testset "9-site DMRG against independent ground spaces" begin
    cases = ((0.0, 11), (0.37, 11))
    @testset "theta=$theta seed=$seed" for (theta, seed) in cases
        check = check_small_system(theta; seed)
        m, result = check.metrics, check.result
        @test m.energy_error_per_site < 1e-8
        @test m.residual_norm < 2e-6
        @test m.ground_space_leakage < 2e-6
        @test m.max_sz_error < 1e-6
        @test m.max_zz_error < 1e-6
        @test m.max_pm_error < 1e-6
        @test m.state_norm_error < 1e-12
        @test m.total_sz ≈ 0.5 atol=1e-12
        @test abs(m.variance) < 1e-10
        @test m.ground_degeneracy == (theta == 0 ? 2 : 1)
        @test length(result.max_truncation_errors) == result.settings.nsweeps
        @test length(result.sweep_energies) == result.settings.nsweeps
        @test all(isfinite, result.max_truncation_errors)
        @test all(>=(0), result.max_truncation_errors)
        @test flux(result.psi) == QN("Sz", 1)
    end
end

@testset "Analytic two-spin truncation controls" begin
    # The field explicitly changes the Hamiltonian. Its Schmidt probabilities
    # are (4/5,1/5); the zero-field singlet has the formerly failing (1/2,1/2).
    cases = ((; name="unequal field control", field=3/8, loss=0.2,
                local_energy=-7/8, retained_energy=-5/8),
             (; name="degenerate singlet", field=0.0, loss=0.5,
                local_energy=-3/4, retained_energy=-1/4))
    @testset "$(case.name)" for case in cases
        (; H, initial) = _heisenberg_chain_control(2; staggered_field=case.field)
        observer = KagomeDMRG._SweepDiagnostics()
        local_energy, truncated = dmrg(H, initial; nsweeps=1, maxdim=1,
            cutoff=0.0, noise=0.0, observer, outputlevel=0)
        @test length(observer.max_truncation_errors) == 1
        @test norm(truncated) ≈ 1 atol=1e-14
        @test only(observer.max_truncation_errors) ≈ case.loss atol=1e-14
        @test local_energy ≈ case.local_energy atol=1e-14
        @test real(inner(truncated', H, truncated)) ≈ case.retained_energy atol=1e-14
    end
end

@testset "Flux-independent longitudinal control" begin
    # This altered Hamiltonian pins a unique product state. It is a readout
    # control, not evidence about the isotropic Heisenberg model.
    lattice = kagome_cylinder(2, 3; Jxy=0, Jz=0)
    sites = spin_sites(lattice)
    labels = [fill("Up", 9); fill("Dn", 9)]
    hz = [fill(1.0, 9); fill(-1.0, 9)]
    baseline = MPS(ComplexF64, sites, labels)
    before = deepcopy(baseline)
    result = run_dmrg(lattice, 0.73; Q=0, sites, psi0=baseline, hz,
                      nsweeps=2, maxdim=2, measure_variance=false)
    @test result.Q == 0 && flux(result.psi) == QN("Sz", 0)
    @test result.energy ≈ -9 atol=1e-12
    @test result.variance === nothing
    @test result.sz ≈ sz_profile(baseline) atol=1e-12
    @test abs(inner(before, baseline)) ≈ 1 atol=1e-12
    @test all(siteind(result.psi, i) == sites[i] for i in 1:18)
    transfer = only(spin_transfer(lattice, sz_profile(baseline), result.sz))
    @test abs(transfer.left) < 1e-12
    @test abs(transfer.right) < 1e-12
    @test abs(transfer.total) < 1e-12
    H0 = twisted_exchange_mpo(sites, lattice, 0; hz)
    @test inner(result.psi', H0, result.psi) ≈ result.energy atol=1e-12

    zero_model = kagome_cylinder(1, 3; Jxy=0, Jz=0)
    zero_sites = spin_sites(zero_model)
    @test norm(full_mpo_matrix(twisted_exchange_mpo(zero_sites, zero_model, 0.2), zero_sites)) == 0
end
