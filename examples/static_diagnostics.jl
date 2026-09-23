#!/usr/bin/env julia
# A bounded current-source fixture, with solve and post-checkpoint measurement
# deliberately run in separate Julia processes. See docs/static_diagnostics.md.
module StaticDiagnosticsExample

using KagomeDMRG
using LinearAlgebra
using SHA
using TOML

const CONFIG_KEYS = Set(["schema_version", "Lx", "Ly", "Q", "seed", "nsweeps", "maxdim"])
const CHECKPOINT_FILES = ("metadata.toml", "state.jls", "checksums.toml")

require(condition, message) = condition || throw(ArgumentError(message))
file_hash(path) = bytes2hex(open(sha256, path))
asdict(value) = Dict(string(k) => v for (k, v) in pairs(value))

function write_toml(path, value)
    open(path, "w") do io
        TOML.print(io, value; sorted=true)
    end
end

function validate_config(cfg)
    require(Set(keys(cfg)) == CONFIG_KEYS, "configuration has missing or unknown keys")
    for key in ("schema_version", "Lx", "Ly", "Q", "seed", "nsweeps")
        require(cfg[key] isa Integer && !(cfg[key] isa Bool), "$key must be an integer")
    end
    require(cfg["schema_version"] == 1, "unsupported configuration schema")
    require(1 <= cfg["Lx"] <= 2 && 3 <= cfg["Ly"] <= 6 &&
        3cfg["Lx"]*cfg["Ly"] <= 18, "example requires a cylinder with at most 18 sites")
    target_sector(3cfg["Lx"]*cfg["Ly"]; Q=cfg["Q"])
    require(0 <= cfg["seed"] <= typemax(Int32), "seed must be a nonnegative Int32")
    require(1 <= cfg["nsweeps"] <= 4, "example permits one to four sweeps")
    dims = cfg["maxdim"]
    require(dims isa AbstractVector && 1 <= length(dims) <= cfg["nsweeps"] &&
        all(d -> d isa Integer && !(d isa Bool) && 1 <= d <= 64, dims),
        "maxdim must be a nonempty schedule with dimensions at most 64")
    return cfg
end

solver_settings(cfg) = (; seed=cfg["seed"], initial_linkdim=4,
    nsweeps=cfg["nsweeps"], maxdim=cfg["maxdim"], cutoff=1e-12, noise=0.0,
    eigsolve_tol=1e-12, eigsolve_krylovdim=20, eigsolve_maxiter=10,
    measure_variance=false)

function require_threads()
    require(Threads.nthreads() == 1, "run the example with --threads=1")
    BLAS.set_num_threads(1)
end

function require_unchanged(pins)
    require(all(isfile(path) && file_hash(path) == hash for (path, hash) in pins),
        "an input, checkpoint, or original solve record changed")
end

"""Create a small NN zero-flux trial checkpoint and a direct measurement reference."""
function solve(output, config_path)
    require_threads()
    cfg = validate_config(TOML.parsefile(config_path))
    output = abspath(output)
    require(!ispath(output), "solve output must be a new directory")
    mkpath(dirname(output))
    mkdir(output)
    write_toml(joinpath(output, "config.toml"), cfg)
    record = Dict{String,Any}("schema_version" => 1, "stage" => "solve",
        "status" => "running", "diagnostic_status" => "not_run",
        "config" => cfg, "driver_sha256" => file_hash(@__FILE__),
        "diagnostic_reference_role" => "direct_measurement_for_process_equivalence",
        "model" => "nearest_neighbor_isotropic_Heisenberg", "theta" => 0.0,
        "gauge" => "seam", "checkpoint_status" => "trial",
        "julia_threads" => 1, "blas_threads" => 1)
    record_path = joinpath(output, "solve.toml")
    write_toml(record_path, record)
    try
        lattice = kagome_cylinder(cfg["Lx"], cfg["Ly"])
        result = run_dmrg(lattice, 0.0; Q=cfg["Q"], solver_settings(cfg)...)
        snapshot = save_checkpoint(output, result;
            baseline=result, theta_path=[0.0], status=:trial)
        # This separate reference is never used to turn a successful solve into
        # a claim that the deferred diagnostic stage or precision checks passed.
        direct = static_diagnostics(result)
        write_toml(joinpath(output, "direct-diagnostics.toml"), direct)
        record["settings"] = asdict(result.settings)
        record["checkpoint"] = relpath(snapshot, output)
        files = vcat(["config.toml", "direct-diagnostics.toml"],
            [joinpath(record["checkpoint"], name) for name in CHECKPOINT_FILES])
        record["file_sha256"] = Dict(name => file_hash(joinpath(output, name)) for name in files)
        require(file_hash(@__FILE__) == record["driver_sha256"], "driver source changed during solve")
        record["status"] = "completed_solve"
        write_toml(record_path, record)
        write(joinpath(output, "solve.sha256"), file_hash(record_path)*"\n")
        return record
    catch err
        record["status"] = "failed"
        record["failure_type"] = string(typeof(err))
        write_toml(record_path, record)
        rethrow()
    end
end

"""Diagnose the pinned trial checkpoint into a fresh directory, preserving its inputs."""
function diagnose(solve_output, output)
    require_threads()
    solve_output, output = abspath(solve_output), abspath(output)
    require(!ispath(output), "diagnostic output must be a new directory")
    record_path = joinpath(solve_output, "solve.toml")
    record_hash = file_hash(record_path)
    require(strip(read(joinpath(solve_output, "solve.sha256"), String)) == record_hash,
        "solve record checksum mismatch")
    record = TOML.parsefile(record_path)
    require(record["schema_version"] == 1 && record["stage"] == "solve" &&
        record["status"] == "completed_solve" && record["diagnostic_status"] == "not_run",
        "expected a completed solve with deferred diagnostics")
    require(record["driver_sha256"] == file_hash(@__FILE__), "driver source differs from solve")
    cfg = validate_config(record["config"])
    require(cfg == TOML.parsefile(joinpath(solve_output, "config.toml")), "solve configuration mismatch")
    require(record["settings"] == asdict(merge(solver_settings(cfg),
        (; initialization="random_fixed_charge"))), "solve settings do not match configuration")
    require(record["model"] == "nearest_neighbor_isotropic_Heisenberg" &&
        record["theta"] == 0.0 && record["gauge"] == "seam" &&
        record["checkpoint_status"] == "trial" &&
        record["julia_threads"] == record["blas_threads"] == 1, "unexpected solve contract")
    parts = splitpath(record["checkpoint"])
    require(length(parts) == 2 && parts[1] == "trial" && startswith(parts[2], "checkpoint-"),
        "checkpoint must belong to the solve output")
    names = Set(vcat(["config.toml", "direct-diagnostics.toml"],
        [joinpath(record["checkpoint"], name) for name in CHECKPOINT_FILES]))
    require(Set(keys(record["file_sha256"])) == names, "incomplete solve artifact pins")
    pins = Dict(joinpath(solve_output, name) => hash for (name, hash) in record["file_sha256"])
    pins[abspath(@__FILE__)] = record["driver_sha256"]
    pins[record_path] = record_hash
    pins[joinpath(solve_output, "solve.sha256")] = file_hash(joinpath(solve_output, "solve.sha256"))
    require_unchanged(pins)
    lattice = kagome_cylinder(cfg["Lx"], cfg["Ly"])
    try
        return diagnose_checkpoint(joinpath(solve_output, record["checkpoint"]), lattice;
            output, Q=cfg["Q"], gauge=:seam, expected_theta=0.0,
            expected_settings=record["settings"], status=:trial)
    finally
        # Repeated diagnostics use separate new directories; even failed runs
        # may never rewrite the solve record or publish a new parent checkpoint.
        require_unchanged(pins)
        if isdir(output)
            write_toml(joinpath(output, "driver.toml"), Dict(
                "schema_version" => 1, "solve_record_sha256" => record_hash,
                "solve_artifacts_unchanged" => true,
                "original_solve_records_modified" => false,
                "driver_sha256" => file_hash(@__FILE__)))
        end
    end
end

function main(args=ARGS)
    require(length(args) == 3 && args[1] in ("solve", "diagnose"),
        "usage: static_diagnostics.jl solve OUTPUT CONFIG | diagnose SOLVE_OUTPUT OUTPUT")
    if args[1] == "solve"
        solve(args[2], args[3])
    else
        diagnose(args[2], args[3])
    end
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    StaticDiagnosticsExample.main()
end
