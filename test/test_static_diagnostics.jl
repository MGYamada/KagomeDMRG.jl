include(joinpath(@__DIR__, "..", "examples", "static_diagnostics.jl"))

@testset "Small static diagnostic driver configuration" begin
    config = TOML.parsefile(joinpath(@__DIR__, "..", "examples", "configs",
                                   "static_diagnostics_small.toml"))
    @test StaticDiagnosticsExample.validate_config(config) === config
    @test_throws ArgumentError StaticDiagnosticsExample.validate_config(
        merge(config, Dict("unknown_option" => true)))
    @test_throws ArgumentError StaticDiagnosticsExample.validate_config(
        merge(config, Dict("maxdim" => [65])))
    @test_throws ArgumentError StaticDiagnosticsExample.validate_config(
        merge(config, Dict("Q" => 0))) # Odd N requires odd integer Q.
    @test_throws ArgumentError StaticDiagnosticsExample.main(String[])
end

# Compare every measurement, including nested Schmidt rows and complex
# correlation components. Timings and process identity are deliberately not
# numerical observables. Tolerances are stricter than the static protocol.
function _static_test_measurements_equal(actual, expected)
    if expected isa AbstractDict
        @test Set(keys(actual)) == Set(keys(expected))
        for key in intersect(keys(actual), keys(expected))
            _static_test_measurements_equal(actual[key], expected[key])
        end
    elseif expected isa AbstractVector{<:Real}
        @test actual ≈ expected atol=1e-10 rtol=0
    elseif expected isa AbstractVector
        @test length(actual) == length(expected)
        for (a, b) in zip(actual, expected)
            _static_test_measurements_equal(a, b)
        end
    elseif expected isa AbstractFloat
        @test actual ≈ expected atol=1e-10 rtol=0
    else
        @test actual == expected
    end
end

_static_test_matrix(rows) = permutedims(hcat(rows...))

@testset "Reusable static diagnostics and immutable post-checkpoint measurement" begin
    lattice = kagome_cylinder(1, 3)
    # One bounded, genuinely optimized fixture; no research-size calculation.
    result = run_dmrg(lattice, 0.0; seed=11, nsweeps=6, maxdim=32,
                      measure_variance=false)
    psi_before = deepcopy(result.psi)
    direct = static_diagnostics(result; chirality=true)
    fixtures = NamedTuple[(; label="optimized", point=result, baseline=result, report=direct)]
    @test direct["schema_version"] == 1
    @test direct["format"] == "KagomeDMRG.static_diagnostics"
    @test direct["reference"] == Dict("provided" => false)
    @test direct["status"] == "completed"
    @test direct["integrity"]["status"] == "passed"
    @test direct["precision"]["status"] == "incomplete"
    @test Set(direct["precision"]["missing"]) ==
          Set(["energy_per_site", "sz_profile", "bond_profile"])
    @test isempty(direct["measurements"]["schmidt"])
    @test result.variance === nothing # Measurement does not rewrite solve results.
    @test all(norm(result.psi[i] - psi_before[i]) == 0 for i in eachindex(result.psi))
    @test direct["measurements"]["energy"] ≈ result.energy atol=1e-10 rtol=0
    @test direct["measurements"]["sz_profile"] ≈ result.sz atol=1e-10 rtol=0
    @test direct["measurements"]["bond_energy"] ≈
          bond_energies(result.psi, lattice) atol=1e-10 rtol=0
    @test direct["measurements"]["variance"] ≈
          real(inner(result.H, result.psi, result.H, result.psi)) - result.energy^2 atol=1e-10 rtol=0
    @test_throws ArgumentError static_diagnostics(result; cuts=[0])
    @test_throws ArgumentError static_diagnostics(result; cuts=[1])

    # A repeated measurement can fill all five comparisons without inventing
    # missing sweep data. This checks status semantics, not new convergence.
    compared = static_diagnostics(result; reference=direct)
    @test compared["precision"]["status"] == "passed"
    @test isempty(compared["precision"]["missing"])
    @test compared["reference"]["provided"]
    reference_bytes = IOBuffer()
    TOML.print(reference_bytes, direct; sorted=true)
    @test compared["reference"]["sha256"] == bytes2hex(sha256(take!(reference_bytes)))
    one_sweep = merge(result, (; settings=merge(result.settings, (; nsweeps=1)),
        sweep_energies=[last(result.sweep_energies)],
        max_truncation_errors=[last(result.max_truncation_errors)]))
    one_report = static_diagnostics(one_sweep; reference=direct)
    @test one_report["precision"]["status"] == "passed"
    @test !one_report["precision"]["within_batch_sweep_change_measured"]
    @test isempty(one_report["precision"]["missing"])
    @test one_report["reference"] == compared["reference"]
    inconsistent = static_diagnostics(merge(result, (; energy=result.energy + 0.1)))
    @test inconsistent["status"] == "completed"
    @test inconsistent["integrity"]["status"] == "failed"
    @test "energy_remeasurement" in inconsistent["integrity"]["failures"]
    @test inconsistent["precision"]["status"] == "not_evaluated"
    @test_throws ArgumentError static_diagnostics(result;
        reference=merge(direct, Dict("theta" => 0.37)))
    for (section, key, value) in (("runtime", "julia", "0.0.0"),
                                  ("provenance", "source_sha256", Dict("changed" => "0"^64)),
                                  ("provenance", "environment_sha256", Dict("changed" => "0"^64)))
        altered = deepcopy(direct)
        altered[section][key] = value
        @test_throws ArgumentError static_diagnostics(result; reference=altered)
    end

    @testset "Independent complex-state reference and unmet precision" begin
        # Synthetic measured state: no claim of optimization. Generic flux,
        # gauge, anisotropy, Q, and field expose hard-coded default conventions.
        model = kagome_cylinder(1, 3; Jxy=0.7, Jz=1.3)
        sites = spin_sites(model)
        psi = initial_mps(sites; Q=-1, seed=43)
        theta, gauge = 0.37, :uniform
        hz = collect(range(-0.2, 0.3; length=9))
        H = twisted_exchange_mpo(sites, model, theta; gauge, hz)
        energy = real(inner(psi', H, psi))
        synthetic = merge(result, (; lattice=model, sites, psi, H, theta, gauge,
            hz, Q=-1, energy, local_energy=energy, sz=sz_profile(psi), variance=nothing,
            sweep_energies=fill(energy, result.settings.nsweeps)))
        report = static_diagnostics(synthetic)
        @test report["integrity"]["status"] == "passed"
        @test report["precision"]["status"] == "unmet"
        @test !report["precision"]["passed"]["variance_per_site"]
        @test !isempty(report["precision"]["missing"])
        @test report["configuration"]["Q"] == -1
        @test report["configuration"]["gauge"] == "uniform"
        @test report["configuration"]["hz"] == hz
        independent = reference_hamiltonian(1, 3, theta; Q=-1, gauge, Jxy=0.7, Jz=1.3)
        vector = sector_amplitudes(psi, sites, independent.basis)
        field_diagonal = [-sum(hz[i] * _reference_sz_bit(s, i) for i in 1:9)
                          for s in independent.basis]
        independent_H = independent.H + Diagonal(field_diagonal)
        independent_energy = real(dot(vector, independent_H * vector))
        second = dot(independent_H * vector, independent_H * vector)
        measured = report["measurements"]
        @test measured["energy"] ≈ independent_energy atol=1e-10 rtol=0
        @test measured["HdaggerH_real"] ≈ real(second) atol=1e-10 rtol=0
        @test measured["HdaggerH_imag"] ≈ imag(second) atol=1e-10 rtol=0
        @test measured["variance"] ≈ real(second)-independent_energy^2 atol=1e-10 rtol=0
        @test measured["variance_per_site"] ≈ measured["variance"]/9 atol=1e-12 rtol=0
        @test measured["sz_profile"] ≈ reference_sz(vector, independent.basis, 9) atol=1e-10 rtol=0
        correlations = measured["correlations"]
        zz = [reference_correlation(vector, independent.basis, 9, i, j) for i in 1:9, j in 1:9]
        pm = [reference_correlation(vector, independent.basis, 9, i, j;
            operators=("S+", "S-")) for i in 1:9, j in 1:9]
        @test maximum(abs, imag.(pm)) > 1e-4
        @test _static_test_matrix(correlations["zz_real"]) ≈ real.(zz) atol=1e-10 rtol=0
        @test _static_test_matrix(correlations["zz_imag"]) ≈ imag.(zz) atol=1e-10 rtol=0
        @test _static_test_matrix(correlations["pm_real"]) ≈ real.(pm) atol=1e-10 rtol=0
        @test _static_test_matrix(correlations["pm_imag"]) ≈ imag.(pm) atol=1e-10 rtol=0
        # Independent gauge angles and correlations, in the public bond order.
        chi = _reference_gauge_angles(1, 3, theta)
        bonds = [b.Jz * real(zz[b.i,b.j]) + b.Jxy *
                 real(cis(b.wy*theta + chi[b.i]-chi[b.j])*pm[b.i,b.j]) for b in model.bonds]
        @test measured["bond_energy"] ≈ bonds atol=1e-10 rtol=0
        @test sum(bonds)-dot(hz, measured["sz_profile"]) ≈ independent_energy atol=1e-10 rtol=0
        # Save the random-state fixture without claiming any optimization. Its
        # zero-flux baseline has the same state/sites/model and an independently
        # calculated zero-flux energy, rather than a reused nonzero-flux value.
        H0 = twisted_exchange_mpo(sites, model, 0.0; gauge, hz)
        independent_H0 = reference_hamiltonian(1, 3, 0.0;
            Q=-1, gauge, Jxy=0.7, Jz=1.3).H + Diagonal(field_diagonal)
        energy0 = real(dot(vector, independent_H0 * vector))
        baseline = merge(synthetic, (; theta=0.0, H=H0, energy=energy0,
            local_energy=energy0, sweep_energies=fill(energy0, result.settings.nsweeps)))
        push!(fixtures, (; label="complex", point=synthetic, baseline, report))
    end

    @testset "Geometric Schmidt cut of an analytic product state" begin
        control = _checkpoint_fixed_field_model(2; Q=0)
        sites = spin_sites(control.lattice)
        psi = MPS(ComplexF64, sites, control.labels)
        H = twisted_exchange_mpo(sites, control.lattice, 0.0; hz=control.hz)
        analytic = merge(result, (; lattice=control.lattice, sites, psi, H,
            hz=control.hz, Q=0, energy=-9.0, local_energy=-9.0,
            sz=sz_profile(psi), variance=nothing,
            sweep_energies=fill(-9.0, result.settings.nsweeps),
            max_truncation_errors=zeros(result.settings.nsweeps)))
        before = deepcopy(psi)
        report = static_diagnostics(analytic)
        @test report["integrity"]["status"] == "passed"
        row = only(report["measurements"]["schmidt"])
        @test row["cut"] == 1
        @test row["bond"] == 9
        @test row["mean_left_sz"] ≈ 4.5 atol=1e-12
        @test row["variance_left_sz"] ≈ 0 atol=1e-12
        @test row["entropy"] ≈ 0 atol=1e-12
        @test report["measurements"]["column_sz"] ≈ [4.5, -4.5] atol=1e-12
        @test all(norm(psi[i] - before[i]) == 0 for i in eachindex(psi))
        @test_throws ArgumentError static_diagnostics(analytic; cuts=[1, 1])
        @test_throws ArgumentError static_diagnostics(analytic; cuts=[2])
        push!(fixtures, (; label="schmidt", point=analytic, baseline=analytic, report))
    end

    mktempdir() do directory
        parent = joinpath(directory, "solve")
        snapshot = save_checkpoint(parent, result; baseline=result,
                                   theta_path=[0.0], status=:trial)
        original = joinpath(parent, "solve.toml")
        write(original, "status = \"completed\"\ndiagnostics_status = \"not_run\"\n")
        before = _checkpoint_test_file_hashes(parent)
        @test !isfile(joinpath(snapshot, "diagnostics.toml"))
        @test TOML.parsefile(original)["diagnostics_status"] == "not_run"

        output = joinpath(directory, "measured")
        saved_cases = NamedTuple[(; label="optimized", snapshot, output, direct,
                                  parent, hashes=before)]
        arguments = ["optimized", snapshot, output]
        for fixture in fixtures[2:end]
            case_parent = joinpath(directory, fixture.label * "-solve")
            path = iszero(fixture.point.theta) ? [0.0] : [0.0, fixture.point.theta]
            case_snapshot = save_checkpoint(case_parent, fixture.point;
                baseline=fixture.baseline, theta_path=path, status=:trial)
            case_output = joinpath(directory, fixture.label * "-measured")
            push!(saved_cases, (; label=fixture.label, snapshot=case_snapshot,
                output=case_output, direct=fixture.report, parent=case_parent,
                hashes=_checkpoint_test_file_hashes(case_parent)))
            append!(arguments, [fixture.label, case_snapshot, case_output])
        end
        worker = joinpath(@__DIR__, "static_diagnostics_worker.jl")
        project = Base.active_project()
        @test project !== nothing
        command = `$(Base.julia_cmd()) --startup-file=no --threads=1 --project=$(dirname(project)) $worker $arguments`
        log = IOBuffer()
        process = run(pipeline(ignorestatus(addenv(command,
            "JULIA_NUM_THREADS" => "1", "OPENBLAS_NUM_THREADS" => "1"));
            stdout=log, stderr=log))
        success(process) || print(String(take!(log)))
        @test success(process)
        if success(process)
            @testset "Fresh-process $(case.label)" for case in saved_cases
                separate = TOML.parsefile(joinpath(case.output, "diagnostics.toml"))
                @test separate["status"] == "completed"
                @test separate["integrity"]["status"] == "passed"
                @test separate["precision"]["status"] == case.direct["precision"]["status"]
                @test separate["configuration"] == case.direct["configuration"]
                @test separate["settings"] == case.direct["settings"]
                _static_test_measurements_equal(separate["measurements"], case.direct["measurements"])
                @test _checkpoint_test_file_hashes(case.parent) == case.hashes
            end
            measured_before = _checkpoint_test_file_hashes(output)
            repeated = diagnose_checkpoint(snapshot, lattice;
                output=joinpath(directory, "remeasured"), chirality=true)
            _static_test_measurements_equal(repeated["measurements"], direct["measurements"])
            @test _checkpoint_test_file_hashes(output) == measured_before
            @test_throws ArgumentError diagnose_checkpoint(snapshot, lattice; output)
            @test _checkpoint_test_file_hashes(output) == measured_before
        end
        @test _checkpoint_test_file_hashes(parent) == before
        @test_throws ArgumentError diagnose_checkpoint(snapshot, lattice;
            output=joinpath(snapshot, "forbidden"))
        @test _checkpoint_test_file_hashes(parent) == before

        failed = joinpath(directory, "incompatible")
        @test_throws ArgumentError diagnose_checkpoint(snapshot, lattice; output=failed, Q=-1)
        failure = TOML.parsefile(joinpath(failed, "diagnostics.toml"))
        @test failure["status"] == "failed"
        @test !haskey(failure, "measurements")
        @test _checkpoint_test_file_hashes(parent) == before
        corrupt = joinpath(directory, "corrupt", "trial", "checkpoint-corrupt")
        mkpath(dirname(corrupt))
        cp(snapshot, corrupt)
        open(joinpath(corrupt, "state.jls"), "a") do io
            write(io, UInt8(0x00))
        end
        corrupt_before = _checkpoint_test_file_hashes(corrupt)
        failed_corrupt = joinpath(directory, "corrupt-diagnostics")
        _checkpoint_test_rejection("checkpoint checksum or byte length mismatch") do
            diagnose_checkpoint(corrupt, lattice; output=failed_corrupt)
        end
        @test TOML.parsefile(joinpath(failed_corrupt, "diagnostics.toml"))["status"] == "failed"
        @test _checkpoint_test_file_hashes(corrupt) == corrupt_before
        @test _checkpoint_test_file_hashes(parent) == before
    end

end
