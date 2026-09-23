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

function _static_workflow_test_config()
    return Dict{String,Any}("schema_version" => 1, "Lx" => 1, "Ly" => 3,
        "Q" => 1, "seed" => 11, "nsweeps" => 1, "maxdim" => [2])
end

function _static_workflow_test_edit_record!(edit, directory)
    path = joinpath(directory, "solve.toml")
    record = TOML.parsefile(path)
    edit(record)
    open(path, "w") do io
        TOML.print(io, record; sorted=true)
    end
    write(joinpath(directory, "solve.sha256"), bytes2hex(sha256(read(path))) * "\n")
    return directory
end

@testset "Shared static configuration contract" begin
    input = _static_workflow_test_config()
    config = static_run_config(input)
    @test config !== input
    @test config["maxdim"] !== input["maxdim"]
    @test config["Jxy"] == config["Jz"] == 1.0
    @test config["hz"] == zeros(9)
    @test config["gauge"] == "seam"
    @test config["initial_linkdim"] == 4
    @test config["cutoff"] == config["eigsolve_tol"] == 1e-12
    @test config["noise"] == 0.0
    @test config["eigsolve_krylovdim"] == 20
    @test config["eigsolve_maxiter"] == 10
    @test static_run_config(config) == config
    input["maxdim"][1] = 7
    @test config["maxdim"] == [2]
    fields = ones(9)
    custom = static_run_config(merge(input, Dict("hz" => fields, "gauge" => "uniform")))
    fields[1] = 3
    @test custom["hz"] == ones(9)

    # Example budgets are not library limits. Validation allocates no MPS.
    large = merge(_static_workflow_test_config(), Dict("Lx" => 3, "Ly" => 9,
        "Q" => 9, "seed" => Int(typemax(Int32)) + 1, "nsweeps" => 5,
        "maxdim" => [128, 256]))
    @test static_run_config(large)["maxdim"] == [128, 256]
    @test_throws ArgumentError static_run_config(merge(input, Dict("unknown_option" => 1)))
    missing = copy(input)
    delete!(missing, "Q")
    @test_throws ArgumentError static_run_config(missing)
    for key in ("schema_version", "Lx", "Ly", "Q", "seed", "nsweeps",
                "initial_linkdim", "cutoff", "noise", "eigsolve_tol",
                "eigsolve_krylovdim", "eigsolve_maxiter", "Jxy", "Jz")
        @test_throws ArgumentError static_run_config(merge(input, Dict(key => true)))
    end
    for (key, value) in (("schema_version", 2), ("Lx", 0), ("Ly", 2),
            ("Q", 0), ("Q", 11), ("seed", -1), ("nsweeps", 0),
            ("initial_linkdim", 0), ("cutoff", -1e-12), ("noise", -1.0),
            ("eigsolve_tol", 0.0), ("eigsolve_krylovdim", 1), ("eigsolve_maxiter", 0),
            ("Jxy", NaN), ("Jz", Inf), ("gauge", "invalid"),
            ("hz", zeros(8)), ("hz", fill(Inf, 9)), ("hz", trues(9)),
            ("maxdim", 2), ("maxdim", Int[]), ("maxdim", [1, 2]),
            ("maxdim", [true]), ("maxdim", [0]), ("maxdim", [2.0]))
        @test_throws ArgumentError static_run_config(merge(input, Dict(key => value)))
    end
end

@testset "Static solve records and deferred diagnosis" begin
    # Zero exchange with a uniform field has E=-Q/2 for every state in this
    # charge sector. One cheap sweep suffices for the workflow contract; this
    # is not a ground-state convergence or measurement-pipeline regression.
    config = merge(_static_workflow_test_config(), Dict("Jxy" => 0.0,
        "Jz" => 0.0, "hz" => ones(9), "initial_linkdim" => 1))
    mktempdir() do directory
        directory = realpath(directory)
        parent = joinpath(directory, "solve")
        events = NamedTuple[]
        solved = run_static(parent, config; progress_callback=event -> push!(events, event))
        record = solved.record
        @test solved.result.energy ≈ -0.5 atol=1e-10 rtol=0
        @test solved.result.variance === nothing
        @test !isempty(events)
        @test record["format"] == "KagomeDMRG.static_run"
        @test record["schema_version"] == 1
        @test record["status"] == "completed_solve"
        @test record["diagnostic_status"] == "not_run"
        @test record["theta"] == 0.0
        @test record["checkpoint_status"] == "trial"
        @test record["config"] == static_run_config(config)
        @test record["settings"]["measure_variance"] === false
        @test record["settings"]["initialization"] == "random_fixed_charge"
        @test record["settings"]["maxdim"] == [2]
        @test record["configuration"]["N"] == 9
        @test record["configuration"]["Q"] == 1
        @test record["configuration"]["hz"] == ones(9)
        @test length(record["configuration"]["sites"]) == 9
        @test !isempty(record["configuration"]["bonds"])
        @test record["runtime"]["julia"] == string(VERSION)
        @test haskey(record["provenance"]["source_sha256"], "src/static_workflow.jl")
        relative_snapshot = record["checkpoint"]
        @test !isabspath(relative_snapshot)
        @test solved.snapshot == joinpath(parent, relative_snapshot)
        expected_files = Set(vcat(["config.toml"],
            [joinpath(relative_snapshot, name)
             for name in ("metadata.toml", "state.jls", "checksums.toml")]))
        @test Set(keys(record["file_sha256"])) == expected_files
        @test all(bytes2hex(sha256(read(joinpath(parent, name)))) == hash
                  for (name, hash) in record["file_sha256"])
        @test !isfile(joinpath(parent, "direct-diagnostics.toml"))
        before = _checkpoint_test_file_hashes(parent)
        loaded = load_static_run(parent)
        @test loaded.config == record["config"]
        @test loaded.record == record
        @test loaded.snapshot == solved.snapshot
        @test loaded.checkpoint.theta == 0.0
        @test loaded.checkpoint.Q == 1
        @test loaded.record_sha256 == strip(read(joinpath(parent, "solve.sha256"), String))
        record["config"]["maxdim"][1] = 33
        @test TOML.parsefile(joinpath(parent, "solve.toml"))["config"]["maxdim"] == [2]
        @test_throws ArgumentError run_static(parent, config)
        @test _checkpoint_test_file_hashes(parent) == before

        @testset "Resealed records still require valid semantics" begin
            alterations = [
                "schema" => r -> (r["schema_version"] = true),
                "format" => r -> (r["format"] = "unrelated"),
                "status" => r -> (r["status"] = "running"),
                "theta" => r -> (r["theta"] = 0.25),
                "path" => r -> (r["checkpoint"] = joinpath("..", relative_snapshot)),
                "config-bool" => r -> (r["config"]["nsweeps"] = true),
                "configuration" => r -> (r["configuration"]["bonds"][1]["wy"] += 1),
                "configuration-bool" => r -> (r["configuration"]["Q"] = true),
                "settings" => r -> (r["settings"]["measure_variance"] = true),
                "settings-bool" => r -> (r["settings"]["initial_linkdim"] = true),
                "runtime" => r -> (r["runtime"]["julia"] = "0.0.0"),
                "provenance" => r -> (r["provenance"]["source_sha256"]["src/static_workflow.jl"] = "0"^64),
                "missing-pin" => r -> delete!(r["file_sha256"], "config.toml"),
                "extra-pin" => r -> (r["file_sha256"]["unexpected.toml"] = "0"^64),
            ]
            for (label, alter) in alterations
                altered = joinpath(directory, label)
                cp(parent, altered)
                _static_workflow_test_edit_record!(alter, altered)
                @test_throws ArgumentError load_static_run(altered)
                if label == "format"
                    rejected_output = joinpath(directory, "malformed-diagnostics")
                    @test_throws ArgumentError diagnose_static_run(altered; output=rejected_output)
                    @test !ispath(rejected_output)
                end
            end
            # Recompute both integrity envelopes: the strict checkpoint loader
            # must still reject a physically incompatible charge claim.
            altered = joinpath(directory, "checkpoint-semantics")
            cp(parent, altered)
            _checkpoint_test_edit_metadata!(joinpath(altered, relative_snapshot)) do metadata
                metadata["state"]["Q"] = -1
            end
            _static_workflow_test_edit_record!(altered) do r
                for name in keys(r["file_sha256"])
                    r["file_sha256"][name] = bytes2hex(sha256(read(joinpath(altered, name))))
                end
            end
            @test_throws ArgumentError load_static_run(altered)
        end

        for name in ("solve.toml", "config.toml", joinpath(relative_snapshot, "state.jls"))
            altered = joinpath(directory, "bytes-" * basename(name))
            cp(parent, altered)
            open(joinpath(altered, name), "a") do io
                write(io, UInt8(0x00))
            end
            @test_throws ArgumentError load_static_run(altered)
        end

        unsealed = joinpath(directory, "unsealed")
        cp(parent, unsealed)
        mv(joinpath(unsealed, "solve.sha256"), joinpath(directory, "removed-seal.sha256"))
        @test_throws ArgumentError load_static_run(unsealed)

        # Matching file bytes cannot authorize an artifact outside its run.
        # Exercise both a single-file alias and a checkpoint-directory alias.
        for (label, name) in (("config-symlink", "config.toml"),
                              ("checkpoint-symlink", relative_snapshot))
            altered = joinpath(directory, label)
            cp(parent, altered)
            local_artifact = joinpath(altered, name)
            external_artifact = joinpath(directory, "external-" * label)
            mv(local_artifact, external_artifact)
            symlink(external_artifact, local_artifact; dir_target=isdir(external_artifact))
            @test_throws ArgumentError load_static_run(altered)
        end

        output = joinpath(directory, "diagnosed")
        report = diagnose_static_run(parent; output)
        @test report["status"] == "completed"
        @test report["integrity"]["status"] == "passed"
        @test report["checkpoint_unchanged"]
        driver = TOML.parsefile(joinpath(output, "driver.toml"))
        @test driver["solve_record_sha256"] == loaded.record_sha256
        @test driver["solve_artifacts_unchanged"]
        @test driver["original_solve_records_modified"] === false
        @test _checkpoint_test_file_hashes(parent) == before
        output_before = _checkpoint_test_file_hashes(output)
        @test_throws ArgumentError diagnose_static_run(parent; output)
        @test _checkpoint_test_file_hashes(output) == output_before
        @test_throws ArgumentError diagnose_static_run(parent; output=joinpath(parent, "forbidden"))
        alias = joinpath(directory, "solve-alias")
        symlink(parent, alias; dir_target=true)
        @test_throws ArgumentError diagnose_static_run(parent; output=joinpath(alias, "forbidden"))
        @test _checkpoint_test_file_hashes(parent) == before

        failed_output = joinpath(directory, "failed-diagnostics")
        @test_throws ArgumentError diagnose_static_run(parent; output=failed_output, cuts=[1])
        failed_report = TOML.parsefile(joinpath(failed_output, "diagnostics.toml"))
        @test failed_report["status"] == "failed"
        @test failed_report["exception_type"] == "ArgumentError"
        @test failed_report["checkpoint_unchanged"]
        failed_driver = TOML.parsefile(joinpath(failed_output, "driver.toml"))
        @test failed_driver["solve_record_sha256"] == loaded.record_sha256
        @test failed_driver["solve_artifacts_unchanged"]
        @test failed_driver["original_solve_records_modified"] === false
        @test _checkpoint_test_file_hashes(parent) == before

        failed = joinpath(directory, "failed")
        private_message = "sentinel exception text must not be persisted"
        @test_throws ErrorException run_static(failed, config;
            progress_callback=event -> error(private_message))
        failure = TOML.parsefile(joinpath(failed, "solve.toml"))
        @test failure["status"] == "failed"
        @test failure["exception_type"] == "ErrorException"
        @test !isempty(failure["active_phase"])
        @test !occursin(private_message, read(joinpath(failed, "solve.toml"), String))
        @test !haskey(failure, "checkpoint")
        @test_throws ArgumentError load_static_run(failed)
        failed_before = _checkpoint_test_file_hashes(failed)
        @test_throws ArgumentError run_static(failed, config)
        @test _checkpoint_test_file_hashes(failed) == failed_before
    end
end
