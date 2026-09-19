using SHA
using TOML

function _continuation_test_policy(; kwargs...)
    thresholds = (; min_overlap=0.0, max_density_change=1.0,
        max_entropy_change=10.0, max_schmidt_change=9.0,
        max_variance=1e-9, max_truncation_error=1e-10,
        max_sweep_energy_change=1e-9, max_cut_spread=1e-9,
        consistency_tol=1e-9)
    return FluxPolicy(; merge(thresholds, (; kwargs...))...)
end

function _continuation_test_hashes(directory)
    return Dict(relpath(joinpath(root, name), directory) =>
        bytes2hex(sha256(read(joinpath(root, name))))
        for (root, _, files) in walkdir(directory) for name in files)
end

@testset "Adaptive sequential flux continuation" begin
    # This deliberately altered Hamiltonian has a unique, flux-independent
    # Q=3 product ground state. Its zero response is known analytically; this
    # control is not a claim about the nearest-neighbor Heisenberg model.
    lattice = kagome_cylinder(3, 3; Jxy=0, Jz=0)
    sites = spin_sites(lattice)
    labels = [fill("Up", 15); fill("Dn", 12)]
    hz = [fill(1.0, 15); fill(-1.0, 12)]
    exact_sz = hz ./ 2
    baseline = run_dmrg(lattice, 0.0; sites,
        psi0=MPS(ComplexF64, sites, labels), hz,
        nsweeps=2, maxdim=2, cutoff=0.0, noise=0.0)
    strict = _continuation_test_policy(min_overlap=1-1e-10,
        max_density_change=1e-10, max_entropy_change=1e-10,
        max_schmidt_change=1e-10)
    @test baseline.energy ≈ -13.5 atol=1e-12
    @test baseline.sz ≈ exact_sz atol=1e-12 rtol=0
    @test flux(baseline.psi) == QN("Sz", 3)

    mktempdir() do directory
        @testset "Analytic zero response, unwrapped cycles, reverse, and restart" begin
            first_leg = continue_flux(lattice, [2pi]; start=baseline,
                output_root=joinpath(directory, "first"), policy=strict,
                initial_step=2pi, min_step=pi/32, hz,
                cuts=[1, 2], diagnostic_bonds=[9, 18], bulk_sites=10:18)
            @test first_leg.status === :completed
            @test first_leg.theta_path == [0.0, 2pi]
            @test length(first_leg.accepted_checkpoints) == 2
            @test length(first_leg.trials) == 1
            @test isfile(first_leg.output_path)
            old_hashes = _continuation_test_hashes(first_leg.last_checkpoint)

            # Restart must retain the measured zero-flux baseline and the
            # earlier unwrapped history, even while returning to theta=0.
            rest = continue_flux(lattice, [4pi, 6pi, 0.0];
                start=first_leg.last_checkpoint,
                output_root=joinpath(directory, "restart"), policy=strict,
                initial_step=2pi, min_step=pi/32, hz,
                cuts=[1, 2], diagnostic_bonds=[9, 18], bulk_sites=10:18)
            @test rest.status === :completed
            @test rest.theta_path == [0.0, 2pi, 4pi, 6pi, 4pi, 2pi, 0.0]
            @test length(rest.trials) == 5
            @test _continuation_test_hashes(first_leg.last_checkpoint) == old_hashes
            @test rest.output_path != first_leg.output_path
            @test isfile(rest.output_path)
            @test TOML.parsefile(rest.output_path)["status"] == "completed"

            for trajectory in (first_leg, rest), trial in trajectory.trials
                @test trial["status"] == "accepted"
                @test isempty(trial["reasons"])
                d = trial["diagnostics"]
                @test d["overlap"] ≈ 1 atol=1e-12
                @test d["max_density_change"] < 1e-12
                @test d["max_entropy_change"] < 1e-12
                @test d["max_schmidt_change"] < 1e-12
                @test abs(d["variance"]) < 1e-10
                @test d["final_truncation_error"] < 1e-12
                @test d["cut_spread"] < 1e-12
                @test d["max_charge_error"] < 1e-12
                @test d["max_schmidt_consistency_error"] < 1e-12
                @test getindex.(d["transfers"], "cut") == [1, 2]
                for transfer in d["transfers"]
                    for key in ("left", "right", "total", "schmidt_left", "schmidt_right")
                        @test abs(transfer[key]) < 1e-12
                    end
                end
                @test getindex.(d["schmidt"], "bond") == [9, 18]
                for (schmidt, charge) in zip(d["schmidt"], (4.5, 6.0))
                    @test schmidt["mean_left_sz"] ≈ charge atol=1e-12
                    @test schmidt["entropy"] ≈ 0 atol=1e-12
                    @test schmidt["probabilities"] ≈ [1.0] atol=1e-12
                end
                # A trial and its accepted publication remain separate files.
                @test !isempty(trial["checkpoint"])
                @test !isempty(trial["accepted_checkpoint"])
                @test trial["checkpoint"] != trial["accepted_checkpoint"]
            end

            for trajectory in (first_leg, rest), path in trajectory.accepted_checkpoints
                checkpoint = load_checkpoint(path, lattice; hz, sites)
                @test checkpoint.diagnostics.energy ≈ -13.5 atol=1e-12
                @test checkpoint.diagnostics.sz ≈ exact_sz atol=1e-12 rtol=0
                @test checkpoint.baseline.sz ≈ exact_sz atol=1e-12 rtol=0
                @test checkpoint.sites == sites
                @test flux(checkpoint.psi) == QN("Sz", 3)
                @test norm(checkpoint.psi) ≈ 1 atol=1e-12
            end
            @test load_checkpoint(rest.last_checkpoint, lattice; hz).theta_path == rest.theta_path
        end

        @testset "Rejected trials reload the accepted state before halving" begin
            calls = NamedTuple[]
            before = Ref{Any}(nothing)
            refinement_root = joinpath(directory, "refinement")
            function reject_once(saved, model, theta; gauge, hz, outputlevel)
                push!(calls, (; from=saved.theta, theta,
                    norm=norm(saved.psi), history=copy(saved.theta_path), sites=copy(saved.sites)))
                if length(calls) == 1
                    accepted_path = only([root for (root, _, files) in walkdir(refinement_root)
                        if basename(dirname(root)) == "accepted" && "metadata.toml" in files])
                    before[] = (; path=accepted_path, hashes=_continuation_test_hashes(accepted_path))
                end
                result = KagomeDMRG._continuation_solve(saved, model, theta;
                    gauge, hz, outputlevel)
                if length(calls) == 1
                    # Corrupt only the in-memory seed after solving. The saved
                    # accepted state must be reloaded rather than reused.
                    saved.psi[1] *= 0
                    return merge(result, (; max_truncation_errors=fill(0.25, result.settings.nsweeps)))
                end
                return result
            end
            refined = continue_flux(lattice, [0.2]; start=baseline,
                output_root=refinement_root, policy=strict,
                initial_step=0.2, min_step=0.025, hz,
                _solver=reject_once)
            @test refined.status === :completed
            @test refined.theta_path == [0.0, 0.1, 0.2]
            @test length(refined.trials) == 3
            @test getindex.(refined.trials, "status") == ["rejected", "accepted", "accepted"]
            @test "final_truncation_error" in refined.trials[1]["reasons"]
            @test !isempty(refined.trials[1]["checkpoint"])
            @test isempty(refined.trials[1]["accepted_checkpoint"])
            @test getproperty.(calls, :from) == [0.0, 0.0, 0.1]
            @test getproperty.(calls, :theta) == [0.2, 0.1, 0.2]
            @test all(c -> isapprox(c.norm, 1; atol=1e-12), calls)
            @test calls[1].history == calls[2].history == [0.0]
            @test all(c -> c.sites == sites, calls)
            @test _continuation_test_hashes(before[].path) == before[].hashes
            trial_path = joinpath(dirname(refined.output_path), refined.trials[1]["checkpoint"])
            rejected = load_checkpoint(trial_path, lattice; hz, status=:trial)
            @test rejected.theta == 0.2
            @test maximum(rejected.diagnostics.max_truncation_errors) == 0.25
            @test_throws ArgumentError load_checkpoint(trial_path, lattice; hz)
        end

        @testset "Numerical failure, minimum step, and trial budget stay unresolved" begin
            calls = NamedTuple[]
            function fail_once(saved, model, theta; gauge, hz, outputlevel)
                push!(calls, (; from=saved.theta, theta, norm=norm(saved.psi)))
                if length(calls) == 1
                    saved.psi[1] *= 0
                    error("injected numerical factorization failure")
                end
                return KagomeDMRG._continuation_solve(saved, model, theta;
                    gauge, hz, outputlevel)
            end
            recovered = continue_flux(lattice, [0.2]; start=baseline,
                output_root=joinpath(directory, "solver-recovery"), policy=strict,
                initial_step=0.2, min_step=0.025, hz, _solver=fail_once)
            @test recovered.status === :completed
            @test recovered.theta_path == [0.0, 0.1, 0.2]
            @test getindex.(recovered.trials, "status") == ["solver_error", "accepted", "accepted"]
            @test recovered.trials[1]["exception_type"] == "ErrorException"
            @test isempty(recovered.trials[1]["checkpoint"])
            @test isempty(recovered.trials[1]["accepted_checkpoint"])
            @test !occursin("injected numerical factorization failure",
                read(recovered.output_path, String))
            @test all(c -> isapprox(c.norm, 1; atol=1e-12), calls)
            @test getproperty.(calls, :from) == [0.0, 0.0, 0.1]

            fail_always(saved, model, theta; kwargs...) = error("injected numerical failure")
            unresolved = continue_flux(lattice, [0.2]; start=baseline,
                output_root=joinpath(directory, "min-step"), policy=strict,
                initial_step=0.2, min_step=0.05, hz, _solver=fail_always)
            @test unresolved.status === :unresolved
            @test unresolved.reason == "min_step"
            @test unresolved.theta_path == [0.0]
            @test length(unresolved.accepted_checkpoints) == 1
            @test getindex.(unresolved.trials, "theta") == [0.2, 0.1, 0.05]
            @test all(t -> t["status"] == "solver_error", unresolved.trials)
            @test TOML.parsefile(unresolved.output_path)["status"] == "unresolved"
            @test load_checkpoint(unresolved.last_checkpoint, lattice; hz).theta == 0.0

            budget = continue_flux(lattice, [0.3]; start=baseline,
                output_root=joinpath(directory, "budget"), policy=strict,
                initial_step=0.1, min_step=0.025, max_trials=1, hz)
            @test budget.status === :unresolved
            @test budget.reason == "max_trials"
            @test budget.theta_path == [0.0, 0.1]
            @test length(budget.trials) == 1
            @test budget.trials[1]["status"] == "accepted"
            @test load_checkpoint(budget.last_checkpoint, lattice; hz).theta == 0.1
        end

        @testset "Invalid trial diagnostics close the journal without promotion" begin
            invalid_root = joinpath(directory, "invalid-trial")
            before = Ref{Any}(nothing)
            function invalid_variance(saved, model, theta; gauge, hz, outputlevel)
                accepted_path = only([root for (root, _, files) in walkdir(invalid_root)
                    if basename(dirname(root)) == "accepted" && "metadata.toml" in files])
                before[] = (; path=accepted_path, hashes=_continuation_test_hashes(accepted_path))
                result = KagomeDMRG._continuation_solve(saved, model, theta;
                    gauge, hz, outputlevel)
                saved.psi[1] *= 0
                return merge(result, (; variance=NaN))
            end
            invalid = continue_flux(lattice, [0.1]; start=baseline,
                output_root=invalid_root, policy=strict,
                initial_step=0.1, min_step=0.025, hz, _solver=invalid_variance)
            @test invalid.status === :unresolved
            @test invalid.reason == "trial_validation_error"
            @test invalid.theta_path == [0.0]
            @test length(invalid.accepted_checkpoints) == 1
            @test invalid.last_checkpoint == before[].path
            @test _continuation_test_hashes(before[].path) == before[].hashes
            @test length(invalid.trials) == 1
            @test invalid.trials[1]["status"] == "invalid_trial"
            @test invalid.trials[1]["reasons"] == ["trial_validation_error"]
            @test invalid.trials[1]["exception_type"] == "ArgumentError"
            @test isempty(invalid.trials[1]["checkpoint"])
            @test isempty(invalid.trials[1]["accepted_checkpoint"])
            @test isempty(invalid.trials[1]["diagnostics"])
            restored = load_checkpoint(invalid.last_checkpoint, lattice; hz, sites)
            @test norm(restored.psi) ≈ 1 atol=1e-12
            @test restored.diagnostics.energy ≈ -13.5 atol=1e-12
            @test restored.diagnostics.sz ≈ exact_sz atol=1e-12 rtol=0
            @test restored.theta == 0.0
            @test restored.theta_path == [0.0]
            journal = TOML.parsefile(invalid.output_path)
            @test journal["status"] == "unresolved"
            @test journal["reason"] == "trial_validation_error"
            @test journal["theta_path"] == [0.0]
            @test length(journal["accepted_checkpoints"]) == 1
            @test journal["trials"][1]["status"] == "invalid_trial"
            @test journal["trials"][1]["exception_type"] == "ArgumentError"
            @test !occursin("checkpoint variance must be finite when measured",
                read(invalid.output_path, String))
        end

        @testset "Policy and protocol input validation" begin
            @test_throws ArgumentError _continuation_test_policy(min_overlap=-0.1)
            @test_throws ArgumentError _continuation_test_policy(min_overlap=1.1)
            for key in (:max_density_change, :max_entropy_change, :max_schmidt_change,
                        :max_variance, :max_truncation_error, :max_sweep_energy_change,
                        :max_cut_spread, :consistency_tol)
                @test_throws ArgumentError _continuation_test_policy(; key => -1.0)
                @test_throws ArgumentError _continuation_test_policy(; key => NaN)
                @test_throws ArgumentError _continuation_test_policy(; key => Inf)
            end
            defaults = (; start=baseline, output_root=joinpath(directory, "invalid"),
                policy=strict, initial_step=0.1, min_step=0.025, hz)
            @test_throws ArgumentError continue_flux(lattice, [NaN]; defaults...)
            @test_throws ArgumentError continue_flux(lattice, [Inf]; defaults...)
            for extra in ((; initial_step=0.0), (; min_step=0.0),
                          (; min_step=0.2), (; max_trials=0),
                          (; cuts=[0]), (; diagnostic_bonds=[0]), (; bulk_sites=[0]),
                          (; cuts=Int[]), (; cuts=[1, 1]),
                          (; bulk_sites=Int[]), (; bulk_sites=[1, 1]),
                          (; gauge=:uniform))
                @test_throws ArgumentError continue_flux(lattice, [0.1]; merge(defaults, extra)...)
            end
            @test_throws ArgumentError continue_flux(lattice, [0.1];
                merge(defaults, (; start=merge(baseline,
                    (; settings=merge(baseline.settings, (; noise=1e-8))))))...)
            @test_throws ArgumentError continue_flux(lattice, [0.1];
                merge(defaults, (; start=merge(baseline, (; variance=nothing,
                    settings=merge(baseline.settings, (; measure_variance=false))))))...)
            @test_throws ArgumentError continue_flux(lattice, [0.1];
                merge(defaults, (; start=merge(baseline, (;
                    settings=merge(baseline.settings, (; nsweeps=1)),
                    sweep_energies=baseline.sweep_energies[1:1],
                    max_truncation_errors=baseline.max_truncation_errors[1:1]))))...)

            # A high-error zero-flux input cannot become accepted merely
            # because it has been supplied as the starting state.
            poor = merge(baseline, (; max_truncation_errors=fill(0.25, 2)))
            initial_failure = continue_flux(lattice, [0.1];
                merge(defaults, (; start=poor,
                    output_root=joinpath(directory, "initial-failure")))...)
            @test initial_failure.status === :unresolved
            @test initial_failure.reason == "initial_diagnostics"
            @test isempty(initial_failure.accepted_checkpoints)
            @test isempty(initial_failure.trials)
            @test initial_failure.last_checkpoint === nothing
            @test TOML.parsefile(initial_failure.output_path)["status"] == "unresolved"

            invalid_initial = continue_flux(lattice, [0.1];
                merge(defaults, (; start=merge(baseline, (; variance=NaN)),
                    output_root=joinpath(directory, "invalid-initial")))...)
            @test invalid_initial.status === :unresolved
            @test invalid_initial.reason == "initial_validation_error"
            @test isempty(invalid_initial.accepted_checkpoints)
            @test isempty(invalid_initial.trials)
            @test invalid_initial.last_checkpoint === nothing
            initial_journal = TOML.parsefile(invalid_initial.output_path)
            @test initial_journal["status"] == "unresolved"
            @test initial_journal["reason"] == "initial_validation_error"
            @test isempty(initial_journal["accepted_checkpoints"])
            @test isempty(initial_journal["trials"])
            @test isempty(initial_journal["last_checkpoint"])
            @test initial_journal["initial"]["status"] == "invalid"
            @test initial_journal["initial"]["exception_type"] == "ArgumentError"
            @test isempty(initial_journal["initial"]["checkpoint"])
            @test isempty(initial_journal["initial"]["diagnostics"])
            @test !occursin("checkpoint variance must be finite when measured",
                read(invalid_initial.output_path, String))

            # Variance can be slightly negative through cancellation, but a
            # macroscopic negative value must not pass a loose magnitude gate.
            negative_variance = continue_flux(lattice, [0.1];
                merge(defaults, (; start=merge(baseline, (; variance=-0.01)),
                    policy=_continuation_test_policy(max_variance=0.1),
                    output_root=joinpath(directory, "negative-variance")))...)
            @test negative_variance.status === :unresolved
            negative_record = TOML.parsefile(negative_variance.output_path)
            @test "negative_variance" in negative_record["initial"]["reasons"]
            @test !("variance" in negative_record["initial"]["reasons"])
            @test negative_record["initial"]["diagnostics"]["negative_variance_tolerance"] < 1e-10
        end
    end
end

@testset "Continuation matches independent nine-site spin-basis ED" begin
    lattice = kagome_cylinder(1, 3)
    baseline = run_dmrg(lattice, 0.0; seed=11)
    mktempdir() do directory
        @test_throws ArgumentError continue_flux(lattice, [0.17]; start=baseline,
            output_root=directory, policy=_continuation_test_policy(),
            initial_step=0.2, min_step=0.0125, diagnostic_bonds=Int[], cuts=Int[])
        trajectory = continue_flux(lattice, [0.17, 0.37]; start=baseline,
            output_root=directory, policy=_continuation_test_policy(),
            initial_step=0.2, min_step=0.0125, diagnostic_bonds=[4], cuts=Int[])
        @test trajectory.status === :completed
        @test trajectory.theta_path ≈ [0.0, 0.17, 0.37] atol=1e-14 rtol=0
        @test length(trajectory.trials) == 2
        for path in trajectory.accepted_checkpoints
            saved = load_checkpoint(path, lattice)
            reference = reference_hamiltonian(1, 3, saved.theta)
            eigensystem = eigen(Hermitian(reference.H))
            v = sector_amplitudes(saved.psi, saved.sites, reference.basis)
            E = saved.diagnostics.energy
            @test abs(E - first(eigensystem.values)) / 9 < 1e-8
            @test norm(reference.H*v-E*v) < 2e-6
            # At zero flux the ground space is twofold degenerate, so use
            # its full projector instead of choosing an arbitrary ED vector.
            indices = findall(e -> abs(e-first(eigensystem.values)) < 1e-9,
                eigensystem.values)
            @test length(indices) == (saved.theta == 0 ? 2 : 1)
            projected = eigensystem.vectors[:, indices] *
                (eigensystem.vectors[:, indices]' * v)
            @test norm(v-projected) < 2e-6
            normalize!(projected)
            @test saved.diagnostics.sz ≈ reference_sz(projected, reference.basis, 9) atol=1e-6 rtol=0
            @test abs(saved.diagnostics.variance) < 1e-9
            @test saved.sites == baseline.sites
            @test saved.baseline.sz == baseline.sz
        end
        @test all(t -> isempty(t["diagnostics"]["transfers"]), trajectory.trials)
    end
end
