@testset "Adaptive sequential flux continuation" begin
    # This deliberately altered Hamiltonian has a unique, flux-independent
    # Q=-3 product ground state. Its zero response is known analytically; this
    # control is not a claim about the nearest-neighbor Heisenberg model.
    baseline = _checkpoint_fixed_field_result(3; Q=-3)
    lattice, sites, hz = baseline.lattice, baseline.sites, baseline.hz
    exact_sz = hz ./ 2
    strict = _checkpoint_test_policy(min_overlap=1-1e-10,
        max_density_change=1e-10, max_entropy_change=1e-10,
        max_schmidt_change=1e-10)
    @test baseline.energy ≈ -13.5 atol=1e-12
    @test baseline.sz ≈ exact_sz atol=1e-12 rtol=0
    @test flux(baseline.psi) == QN("Sz", -3)

    mktempdir() do directory
        @testset "Analytic zero response, unwrapped cycles, reverse, and restart" begin
            first_leg = continue_flux(lattice, [2pi]; start=baseline,
                output_root=joinpath(directory, "first"), policy=strict,
                initial_step=2pi, min_step=pi/32, hz, Q=-3,
                cuts=[1, 2], diagnostic_bonds=[9, 18], bulk_sites=10:18)
            @test first_leg.status === :completed
            @test first_leg.theta_path == [0.0, 2pi]
            @test length(first_leg.accepted_checkpoints) == 2
            @test length(first_leg.trials) == 1
            @test isfile(first_leg.output_path)
            old_hashes = _checkpoint_test_file_hashes(first_leg.last_checkpoint)

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
            @test _checkpoint_test_file_hashes(first_leg.last_checkpoint) == old_hashes
            @test rest.output_path != first_leg.output_path
            @test isfile(rest.output_path)
            @test TOML.parsefile(rest.output_path)["status"] == "completed"

            # Inspect one representative trial, rather than repeating every
            # diagnostic for all six steps of the same exact product trajectory.
            trial = last(rest.trials)
            d = trial["diagnostics"]
            @test trial["status"] == "accepted"
            @test isempty(trial["reasons"])
            @test d["overlap"] ≈ 1 atol=1e-12
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
            for (schmidt, charge) in zip(d["schmidt"], (4.5, 3.0))
                @test schmidt["mean_left_sz"] ≈ charge atol=1e-12
                @test schmidt["probabilities"] ≈ [1.0] atol=1e-12
            end
            @test trial["checkpoint"] != trial["accepted_checkpoint"]

            # The baseline above and the restored final state cover both ends.
            checkpoint = load_checkpoint(rest.last_checkpoint, lattice; hz, sites)
            @test checkpoint.Q == -3
            @test checkpoint.diagnostics.energy ≈ -13.5 atol=1e-12
            @test checkpoint.diagnostics.sz ≈ exact_sz atol=1e-12 rtol=0
            @test checkpoint.baseline.sz ≈ exact_sz atol=1e-12 rtol=0
            @test checkpoint.theta_path == rest.theta_path

        end

        @testset "Rejected trials reload the accepted state before halving" begin
            # The analytic control solver isolates rollback and step control;
            # the full forward/reverse trajectory above uses actual DMRG.
            calls = NamedTuple[]
            before = Ref{Any}(nothing)
            refinement_root = joinpath(directory, "refinement")
            function reject_once(saved, model, theta; gauge, hz, outputlevel)
                push!(calls, (; from=saved.theta, theta,
                    norm=norm(saved.psi), history=copy(saved.theta_path), sites=copy(saved.sites)))
                if length(calls) == 1
                    accepted_path = only([root for (root, _, files) in walkdir(refinement_root)
                        if basename(dirname(root)) == "accepted" && "metadata.toml" in files])
                    before[] = (; path=accepted_path, hashes=_checkpoint_test_file_hashes(accepted_path))
                end
                result = _checkpoint_fixed_field_step(saved, model, theta;
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
            @test _checkpoint_test_file_hashes(before[].path) == before[].hashes
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
                return _checkpoint_fixed_field_step(saved, model, theta;
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
                initial_step=0.1, min_step=0.025, max_trials=1, hz,
                _solver=_checkpoint_fixed_field_step)
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
                before[] = (; path=accepted_path, hashes=_checkpoint_test_file_hashes(accepted_path))
                result = _checkpoint_fixed_field_step(saved, model, theta;
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
            @test _checkpoint_test_file_hashes(before[].path) == before[].hashes
            @test length(invalid.trials) == 1
            @test invalid.trials[1]["status"] == "invalid_trial"
            @test invalid.trials[1]["reasons"] == ["trial_validation_error"]
            @test isempty(invalid.trials[1]["checkpoint"])
            @test isempty(invalid.trials[1]["accepted_checkpoint"])
            @test isempty(invalid.trials[1]["diagnostics"])
            restored = load_checkpoint(invalid.last_checkpoint, lattice; hz, sites)
            @test norm(restored.psi) ≈ 1 atol=1e-12
            @test restored.theta_path == [0.0]
            journal = TOML.parsefile(invalid.output_path)
            @test journal["status"] == "unresolved"
            @test journal["reason"] == "trial_validation_error"
            @test journal["trials"][1]["exception_type"] == "ArgumentError"
            @test !occursin("checkpoint variance must be finite when measured",
                read(invalid.output_path, String))
        end

        @testset "Policy and protocol input validation" begin
            @testset "Invalid overlap $value" for value in (-0.1, 1.1)
                @test_throws ArgumentError _checkpoint_test_policy(min_overlap=value)
            end
            # One invalid representative per field exercises the shared
            # finite/nonnegative validator without a field × value product.
            @testset "Invalid $key = $value" for (key, value) in (
                    (:max_density_change, -1.0), (:max_entropy_change, NaN),
                    (:max_schmidt_change, Inf), (:max_variance, -1.0),
                    (:max_truncation_error, NaN), (:max_sweep_energy_change, Inf),
                    (:max_cut_spread, -1.0), (:consistency_tol, NaN))
                @test_throws ArgumentError _checkpoint_test_policy(; key => value)
            end
            defaults = (; start=baseline, output_root=joinpath(directory, "invalid"),
                policy=strict, initial_step=0.1, min_step=0.025, hz)
            @testset "Nonfinite target $theta" for theta in (NaN, Inf)
                @test_throws ArgumentError continue_flux(lattice, [theta]; defaults...)
            end
            @testset "Invalid protocol: $label" for (label, extra) in (
                ("zero initial step", (; initial_step=0.0)),
                ("zero minimum step", (; min_step=0.0)),
                ("minimum above initial", (; min_step=0.2)),
                ("zero trial budget", (; max_trials=0)),
                ("different charge", (; Q=3)),
                ("noninteger charge", (; Q=-3.0)),
                ("cut index", (; cuts=[0])), ("bond index", (; diagnostic_bonds=[0])),
                ("bulk index", (; bulk_sites=[0])), ("no cuts", (; cuts=Int[])),
                ("repeated cuts", (; cuts=[1, 1])), ("no bulk sites", (; bulk_sites=Int[])),
                ("repeated bulk sites", (; bulk_sites=[1, 1])), ("gauge", (; gauge=:uniform)))
                @test_throws ArgumentError continue_flux(lattice, [0.1]; merge(defaults, extra)...)
            end
            @testset "Invalid starting solver: $label" for (label, changed) in (
                ("noise", (; settings=merge(baseline.settings, (; noise=1e-8)))),
                ("unmeasured variance", (; variance=nothing,
                    settings=merge(baseline.settings, (; measure_variance=false)))),
                ("one sweep", (; settings=merge(baseline.settings, (; nsweeps=1)),
                    sweep_energies=baseline.sweep_energies[1:1],
                    max_truncation_errors=baseline.max_truncation_errors[1:1])))
                @test_throws ArgumentError continue_flux(lattice, [0.1];
                    merge(defaults, (; start=merge(baseline, changed)))...)
            end

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
                    policy=_checkpoint_test_policy(max_variance=0.1),
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
            output_root=directory, policy=_checkpoint_test_policy(),
            initial_step=0.2, min_step=0.0125, diagnostic_bonds=Int[], cuts=Int[])
        trajectory = continue_flux(lattice, [0.17, 0.37]; start=baseline,
            output_root=directory, policy=_checkpoint_test_policy(),
            initial_step=0.2, min_step=0.0125, diagnostic_bonds=[4], cuts=Int[])
        @test trajectory.status === :completed
        @test trajectory.theta_path ≈ [0.0, 0.17, 0.37] atol=1e-14 rtol=0
        @test length(trajectory.trials) == 2
        # Zero flux checks the degenerate projector; the generic endpoint
        # checks the unique ground state. The intermediate point is redundant.
        for path in (first(trajectory.accepted_checkpoints), trajectory.last_checkpoint)
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
