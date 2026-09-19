@testset "Atomic local checkpoints and validated restart" begin
    lattice = kagome_cylinder(1, 3)
    # Reuse the existing solves in a nondefault sector so restart must inherit Q.
    baseline = run_dmrg(lattice, 0.0; seed=11, Q=-1)
    point = run_dmrg(lattice, 0.37; sites=baseline.sites, psi0=baseline.psi, seed=11, Q=-1)
    @test point.lattice.sites !== lattice.sites
    @test point.lattice.bonds !== lattice.bonds
    @test point.hz == zeros(9)

    mktempdir() do directory
        accepted = save_checkpoint(directory, point; baseline,
                                   theta_path=[0.0, 0.37], status=:accepted)
        @test basename(dirname(accepted)) == "accepted"
        @test all(isfile(joinpath(accepted, name)) for name in
                  ("metadata.toml", "state.jls", "checksums.toml"))
        before = _checkpoint_test_file_hashes(accepted)
        loaded = load_checkpoint(accepted, lattice; sites=point.sites,
                                 expected_theta=0.37, expected_settings=point.settings)
        @test loaded.theta == 0.37
        @test loaded.theta_path == [0.0, 0.37]
        @test loaded.gauge === :seam
        @test loaded.Q == -1
        @test loaded.hz == zeros(9)
        @test loaded.settings == point.settings
        @test loaded.execution_identity == point.execution_identity
        @test loaded.baseline.execution_identity == baseline.execution_identity
        @test loaded.sites == point.sites
        @test all(siteind(loaded.psi, i) == point.sites[i] for i in 1:9)
        @test flux(loaded.psi) == QN("Sz", -1)
        @test norm(loaded.psi) ≈ 1 atol=2e-12
        @test abs(inner(loaded.psi, point.psi)) ≈ 1 atol=2e-12
        @test sz_profile(loaded.psi) ≈ point.sz atol=2e-12 rtol=0
        @test loaded.baseline.theta == 0.0
        @test loaded.baseline.sz == baseline.sz
        @test loaded.baseline.sites == baseline.sites
        @test loaded.baseline.Q == -1
        @test abs(inner(loaded.baseline.psi, baseline.psi)) ≈ 1 atol=2e-12
        second = save_checkpoint(directory, point; baseline,
                                 theta_path=[0.0, 0.37], status=:accepted)
        @test second != accepted
        @test _checkpoint_test_file_hashes(accepted) == before

        # Change H while retaining an explicitly unwrapped theta. Reusing the
        # old MPO or reducing theta modulo 2pi cannot satisfy both checks.
        target = 2pi + 0.71
        solver = (; (k => v for (k, v) in pairs(point.settings)
                      if k != :initialization)...)
        direct = run_dmrg(lattice, target; sites=point.sites, psi0=point.psi, Q=-1, solver...)
        resumed = resume_dmrg(accepted, lattice, target;
                              expected_settings=point.settings)
        @test resumed.theta == target
        @test resumed.settings.initialization == "provided_mps"
        @test resumed.execution_identity == point.execution_identity
        @test resumed.sites == direct.sites
        @test resumed.energy ≈ direct.energy atol=1e-10 rtol=0
        @test resumed.sz ≈ direct.sz atol=1e-7 rtol=0
        @test abs(inner(resumed.psi, direct.psi)) ≈ 1 atol=1e-10
        reference = reference_hamiltonian(1, 3, target; Q=-1)
        v = sector_amplitudes(resumed.psi, resumed.sites, reference.basis)
        @test norm(reference.H * v - resumed.energy * v) < 2e-6
        @test abs(resumed.energy - eigmin(Hermitian(reference.H))) / 9 < 1e-8

        # Trial writes never replace an accepted snapshot. A reloaded baseline
        # can seed a fresh checkpoint without rebuilding a zero-flux state.
        trial = save_checkpoint(directory, resumed; baseline=loaded.baseline,
                                theta_path=[0.0, 0.37, target])
        @test basename(dirname(trial)) == "trial"
        @test _checkpoint_test_file_hashes(accepted) == before
        @test_throws ArgumentError load_checkpoint(trial, lattice)
        @test_throws ArgumentError resume_dmrg(trial, lattice, target + 0.1)
        trial_loaded = load_checkpoint(trial, lattice; status=:trial, Q=-1)
        @test trial_loaded.theta_path == [0.0, 0.37, target]
        @test trial_loaded.baseline.sz == baseline.sz
        reverse_path = save_checkpoint(directory, point; baseline,
            theta_path=[0.0, target, 0.37])
        @test load_checkpoint(reverse_path, lattice; status=:trial).theta_path ==
              [0.0, target, 0.37]

        @testset "Reject incompatible requests and invalid snapshots" begin
            shifted = deepcopy(lattice)
            s = shifted.sites[1]
            shifted.sites[1] = KagomeSite(s.index, s.x, s.y, s.sublattice,
                                         (s.position[1] + 0.01, s.position[2]))
            for (label, model, options) in (
                ("geometry size", kagome_cylinder(2, 3), (;)),
                ("coupling", kagome_cylinder(1, 3; Jxy=0.9), (;)),
                ("site position", shifted, (;)),
                ("gauge", lattice, (; gauge=:uniform)),
                ("charge", lattice, (; Q=1)),
                ("noninteger charge", lattice, (; Q=-1.0)),
                ("field", lattice, (; hz=ones(9))),
                ("theta", lattice, (; expected_theta=0.4)),
                ("site identities", lattice, (; sites=spin_sites(lattice))),
                ("settings", lattice, (; expected_settings=merge(point.settings, (; cutoff=1e-9)))),
                ("status", lattice, (; status=:trial)))
                @testset "Incompatible $label" begin
                    @test_throws ArgumentError load_checkpoint(accepted, model; options...)
                end
            end
            @test_throws ArgumentError save_checkpoint(directory, point; baseline,
                theta_path=[0.0, 0.37], status=:invalid)
            @test_throws ArgumentError resume_dmrg(accepted, lattice, target; Q=1)
            missing_identity = (; (k=>v for (k,v) in pairs(point)
                                   if k != :execution_identity)...)
            _checkpoint_test_rejection("missing its execution identity") do
                save_checkpoint(directory, missing_identity; baseline,
                    theta_path=[0.0, 0.37])
            end
            different_identity = KagomeDMRG._ExecutionIdentity("different.toml",
                point.execution_identity.source_sha256, point.execution_identity.runtime)
            _checkpoint_test_rejection("execution identity does not match") do
                save_checkpoint(directory, merge(point, (; execution_identity=different_identity));
                    baseline, theta_path=[0.0, 0.37])
            end
            _checkpoint_test_rejection("execution identity does not match") do
                save_checkpoint(directory, point;
                    baseline=merge(baseline, (; execution_identity=different_identity)),
                    theta_path=[0.0, 0.37])
            end
            for (label, theta_path) in (
                ("missing zero", [0.37]), ("wrapped endpoint", [0.0, 0.37 + 2pi]),
                ("nonfinite", [0.0, NaN, 0.37]))
                @testset "Invalid path: $label" begin
                    @test_throws ArgumentError save_checkpoint(directory, point; baseline, theta_path)
                end
            end
            @test_throws ArgumentError save_checkpoint(directory, point;
                baseline=merge(baseline, (; theta=0.1)), theta_path=[0.0, 0.37])
            @test_throws ArgumentError save_checkpoint(directory, point;
                baseline=merge(baseline, (; Q=1)), theta_path=[0.0, 0.37])
            bad_norm = deepcopy(point.psi)
            bad_norm[1] *= 2
            for (label, changed) in (
                ("charge", (; Q=3)), ("profile", (; sz=point.sz .+ 0.01)),
                ("site identities", (; sites=spin_sites(lattice))), ("norm", (; psi=bad_norm)))
                @testset "Invalid result: $label" begin
                    @test_throws ArgumentError save_checkpoint(directory,
                        merge(point, changed); baseline, theta_path=[0.0, 0.37])
                end
            end
            _checkpoint_test_rejection("checkpoint energy does not match its MPS and Hamiltonian") do
                save_checkpoint(directory, merge(point, (; energy=point.energy + 0.1));
                    baseline, theta_path=[0.0, 0.37])
            end
            phase_changed = deepcopy(point.psi)
            phase_changed[1] = noprime(2 * op("Sz", point.sites[1]) * phase_changed[1])
            @test sz_profile(phase_changed) ≈ point.sz atol=2e-12 rtol=0
            @test norm(phase_changed) ≈ 1 atol=2e-12
            @test abs(real(inner(phase_changed', point.H, phase_changed))-point.energy) > 1e-3
            _checkpoint_test_rejection("checkpoint energy does not match its MPS and Hamiltonian") do
                save_checkpoint(directory, merge(point, (; psi=phase_changed));
                    baseline, theta_path=[0.0, 0.37])
            end
            profile_only_baseline = (; (k=>v for (k,v) in pairs(loaded.baseline) if k != :psi)...)
            _checkpoint_test_rejection("baseline must include its measured zero-flux MPS") do
                save_checkpoint(directory, point; baseline=profile_only_baseline,
                    theta_path=[0.0, 0.37])
            end
            @test _checkpoint_test_file_hashes(accepted) == before
        end

        @testset "Integrity envelope and semantic metadata validation" begin
            @testset "Truncated $filename" for filename in ("metadata.toml", "state.jls")
                copied = joinpath(directory, "accepted", "checkpoint-truncated-" * filename)
                cp(accepted, copied)
                open(joinpath(copied, filename), "w") do io
                    write(io, UInt8[0x00, 0x01])
                end
                _checkpoint_test_rejection("checkpoint checksum or byte length mismatch") do
                    load_checkpoint(copied, lattice)
                end
            end
            copied = joinpath(directory, "accepted", "checkpoint-missing-checksums")
            cp(accepted, copied)
            rm(joinpath(copied, "checksums.toml"))
            _checkpoint_test_rejection("invalid checkpoint") do
                load_checkpoint(copied, lattice)
            end

            @testset "Forged metadata: $label" for (label, edit, message) in (
                ("schema", m -> (m["schema_version"] = 999), "unsupported checkpoint schema"),
                ("format", m -> (m["format"] = "unrecognized"), "unsupported checkpoint schema"),
                ("runtime", m -> (m["runtime"]["julia"] = "0.0.0"), "checkpoint runtime mismatch"),
                ("source", m -> (m["provenance"]["source_sha256"]["src/model.jl"] = "0"^64),
                 "checkpoint source or dependency manifest mismatch"),
                ("charge", m -> (m["state"]["Q"] = 3), "checkpoint charge metadata mismatch"),
                ("configuration-charge-type", m -> (m["configuration"]["Q"] = -1.0),
                    "checkpoint charge must be an integer"),
                ("bond-family", m -> (m["configuration"]["bonds"][1]["family"] = "J3"),
                    "checkpoint model or geometry mismatch"),
                ("state-charge-type", m -> (m["state"]["Q"] = -1.0),
                 "checkpoint charge must be an integer"),
                ("path", m -> (m["theta_path"] = [0.0, 0.99]), "checkpoint flux path must"),
                ("energy", m -> (m["state"]["energy"] += 0.1),
                 "checkpoint energy does not match its MPS and Hamiltonian"),
                ("profile", m -> (m["state"]["sz"][1:2] .+= [0.001, -0.001]),
                 "checkpoint stored profile differs from MPS"),
                ("baseline-profile", m -> (m["baseline"]["sz"][1:2] .+= [0.001, -0.001]),
                 "baseline profile does not match the zero-flux MPS"))
                copied = joinpath(directory, "accepted", "checkpoint-wrong-" * label)
                cp(accepted, copied)
                _checkpoint_test_edit_metadata!(edit, copied)
                _checkpoint_test_rejection(message) do
                    load_checkpoint(copied, lattice)
                end
            end
            copied = joinpath(directory, "accepted", "checkpoint-checksum-schema")
            cp(accepted, copied)
            checksum_path = joinpath(copied, "checksums.toml")
            checksums = TOML.parsefile(checksum_path)
            checksums["schema_version"] = 999
            open(checksum_path, "w") do io
                TOML.print(io, checksums; sorted=true)
            end
            _checkpoint_test_rejection("unsupported checksum schema") do
                load_checkpoint(copied, lattice)
            end
            @testset "Invalid serialized state: $label" for label in ("norm", "charge", "phase")
                copied = joinpath(directory, "accepted", "checkpoint-state-" * label)
                cp(accepted, copied)
                psi = if label == "norm"
                    state = deepcopy(point.psi)
                    state[1] *= 2
                    state
                elseif label == "charge"
                    MPS(ComplexF64, point.sites, fill("Up", 9))
                else
                    state = deepcopy(point.psi)
                    state[1] = noprime(2 * op("Sz", point.sites[1]) * state[1])
                    state
                end
                payload = open(deserialize, joinpath(copied, "state.jls"))
                open(joinpath(copied, "state.jls"), "w") do io
                    serialize(io, merge(payload, (; psi)))
                end
                _checkpoint_test_rehash!(copied, "state.jls")
                expected = label == "norm" ? "checkpoint state must be normalized" :
                           label == "charge" ? "checkpoint total charge mismatch" :
                           "checkpoint energy does not match its MPS and Hamiltonian"
                _checkpoint_test_rejection(expected) do
                    load_checkpoint(copied, lattice)
                end
            end
            @test _checkpoint_test_file_hashes(accepted) == before
        end
    end
end

@testset "Zero-charge checkpoint on a lattice outside the default 1/9 sizes" begin
    # This exact product state needs no optimizer. Its model/profile/energy are
    # still checked by save/load, including the even N=12, Q=0 sector contract.
    lattice = kagome_cylinder(1, 4; Jxy=0.0, Jz=0.0)
    sites = spin_sites(lattice)
    hz = [fill(1.0, 6); fill(-1.0, 6)]
    psi = MPS(ComplexF64, sites, [fill("Up", 6); fill("Dn", 6)])
    settings = (; seed=0, initial_linkdim=1, nsweeps=2, maxdim=[1], cutoff=0.0,
        noise=0.0, eigsolve_tol=1e-12, eigsolve_krylovdim=3, eigsolve_maxiter=1,
        measure_variance=true, initialization="provided_mps")
    point = (; lattice, sites, hz, psi, settings, theta=0.0, gauge=:seam, Q=0,
        execution_identity=KagomeDMRG._execution_identity(), energy=-6.0,
        local_energy=-6.0, variance=0.0, sz=hz ./ 2,
        sweep_energies=fill(-6.0, 2), max_truncation_errors=zeros(2))
    mktempdir() do directory
        path = save_checkpoint(directory, point; baseline=point, theta_path=[0.0], status=:accepted)
        @test load_checkpoint(path, lattice; hz).Q == 0
        @test_throws ArgumentError save_checkpoint(directory, merge(point, (; Q=1));
            baseline=point, theta_path=[0.0])
    end
end
