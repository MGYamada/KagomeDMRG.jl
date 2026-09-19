using SHA
using Serialization
using TOML

function _checkpoint_test_file_hashes(directory)
    return Dict(relpath(joinpath(root, name), directory) =>
                bytes2hex(sha256(read(joinpath(root, name))))
                for (root, _, files) in walkdir(directory) for name in files)
end

# Recompute the public integrity envelope after an intentional data change.
# This distinguishes semantic validation from detection of accidental corruption.
function _checkpoint_test_rehash!(directory, filename)
    path = joinpath(directory, filename)
    checksum_path = joinpath(directory, "checksums.toml")
    checksums = TOML.parsefile(checksum_path)
    checksums["files"][filename]["sha256"] = bytes2hex(sha256(read(path)))
    checksums["files"][filename]["bytes"] = filesize(path)
    open(checksum_path, "w") do io
        TOML.print(io, checksums; sorted=true)
    end
    return directory
end

function _checkpoint_test_edit_metadata!(edit, directory)
    path = joinpath(directory, "metadata.toml")
    metadata = TOML.parsefile(path)
    edit(metadata)
    open(path, "w") do io
        TOML.print(io, metadata; sorted=true)
    end
    return _checkpoint_test_rehash!(directory, "metadata.toml")
end

function _checkpoint_test_rejection(action, expected_message)
    exception = try
        action()
        nothing
    catch error
        error
    end
    @test exception isa ArgumentError
    if exception isa Exception
        @test occursin(expected_message, sprint(showerror, exception))
    end
end

@testset "Atomic local checkpoints and validated restart" begin
    lattice = kagome_cylinder(1, 3)
    baseline = run_dmrg(lattice, 0.0; seed=11)
    point = run_dmrg(lattice, 0.37; sites=baseline.sites, psi0=baseline.psi, seed=11)
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
        @test loaded.Q == 1
        @test loaded.hz == zeros(9)
        @test loaded.settings == point.settings
        @test loaded.sites == point.sites
        @test all(siteind(loaded.psi, i) == point.sites[i] for i in 1:9)
        @test flux(loaded.psi) == QN("Sz", 1)
        @test norm(loaded.psi) ≈ 1 atol=2e-12
        @test abs(inner(loaded.psi, point.psi)) ≈ 1 atol=2e-12
        @test sz_profile(loaded.psi) ≈ point.sz atol=2e-12 rtol=0
        @test loaded.baseline.theta == 0.0
        @test loaded.baseline.sz == baseline.sz
        @test loaded.baseline.sites == baseline.sites
        @test loaded.baseline.Q == 1
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
        direct = run_dmrg(lattice, target; sites=point.sites, psi0=point.psi, solver...)
        resumed = resume_dmrg(accepted, lattice, target;
                              expected_settings=point.settings)
        @test resumed.theta == target
        @test resumed.settings.initialization == "provided_mps"
        @test resumed.sites == direct.sites
        @test resumed.energy ≈ direct.energy atol=1e-10 rtol=0
        @test resumed.sz ≈ direct.sz atol=1e-7 rtol=0
        @test abs(inner(resumed.psi, direct.psi)) ≈ 1 atol=1e-10
        reference = reference_hamiltonian(1, 3, target)
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
        trial_loaded = load_checkpoint(trial, lattice; status=:trial)
        @test trial_loaded.theta_path == [0.0, 0.37, target]
        @test trial_loaded.baseline.sz == baseline.sz
        reverse_path = save_checkpoint(directory, point; baseline,
            theta_path=[0.0, target, 0.37])
        @test load_checkpoint(reverse_path, lattice; status=:trial).theta_path ==
              [0.0, target, 0.37]

        @testset "Reject incompatible requests and invalid snapshots" begin
            @test_throws ArgumentError load_checkpoint(accepted, kagome_cylinder(2, 3))
            @test_throws ArgumentError load_checkpoint(accepted,
                                                       kagome_cylinder(1, 3; Jxy=0.9))
            shifted = deepcopy(lattice)
            s = shifted.sites[1]
            shifted.sites[1] = KagomeSite(s.index, s.x, s.y, s.sublattice,
                                         (s.position[1] + 0.01, s.position[2]))
            @test_throws ArgumentError load_checkpoint(accepted, shifted)
            @test_throws ArgumentError load_checkpoint(accepted, lattice; gauge=:uniform)
            @test_throws ArgumentError load_checkpoint(accepted, lattice; hz=ones(9))
            @test_throws ArgumentError load_checkpoint(accepted, lattice; expected_theta=0.4)
            @test_throws ArgumentError load_checkpoint(accepted, lattice; sites=spin_sites(lattice))
            @test_throws ArgumentError load_checkpoint(accepted, lattice;
                expected_settings=merge(point.settings, (; cutoff=1e-9)))
            @test_throws ArgumentError load_checkpoint(accepted, lattice; status=:trial)
            @test_throws ArgumentError save_checkpoint(directory, point; baseline,
                theta_path=[0.0, 0.37], status=:invalid)
            @test_throws ArgumentError save_checkpoint(directory, point; baseline,
                theta_path=[0.37])
            @test_throws ArgumentError save_checkpoint(directory, point; baseline,
                theta_path=[0.0, 0.37 + 2pi])
            @test_throws ArgumentError save_checkpoint(directory, point; baseline,
                theta_path=[0.0, NaN, 0.37])
            @test_throws ArgumentError save_checkpoint(directory, point;
                baseline=merge(baseline, (; theta=0.1)), theta_path=[0.0, 0.37])
            @test_throws ArgumentError save_checkpoint(directory,
                merge(point, (; Q=3)); baseline, theta_path=[0.0, 0.37])
            @test_throws ArgumentError save_checkpoint(directory,
                merge(point, (; sz=point.sz .+ 0.01)); baseline, theta_path=[0.0, 0.37])
            @test_throws ArgumentError save_checkpoint(directory,
                merge(point, (; sites=spin_sites(lattice))); baseline,
                theta_path=[0.0, 0.37])
            bad_norm = deepcopy(point.psi)
            bad_norm[1] *= 2
            @test_throws ArgumentError save_checkpoint(directory,
                merge(point, (; psi=bad_norm)); baseline, theta_path=[0.0, 0.37])
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
            for filename in ("metadata.toml", "state.jls")
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

            for (label, edit, message) in (
                ("schema", m -> (m["schema_version"] = 999), "unsupported checkpoint schema"),
                ("format", m -> (m["format"] = "unrecognized"), "unsupported checkpoint schema"),
                ("runtime", m -> (m["runtime"]["julia"] = "0.0.0"), "checkpoint runtime mismatch"),
                ("source", m -> (m["provenance"]["source_sha256"]["src/model.jl"] = "0"^64),
                 "checkpoint source or dependency manifest mismatch"),
                ("charge", m -> (m["state"]["Q"] = 3), "checkpoint charge metadata mismatch"),
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
            for label in ("norm", "charge", "phase")
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
