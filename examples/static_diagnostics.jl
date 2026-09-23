#!/usr/bin/env julia
# Bounded example of the shared static-run contract. See docs/static_workflow.md.
module StaticDiagnosticsExample

using KagomeDMRG
using LinearAlgebra
using TOML

require(condition, message) = condition || throw(ArgumentError(message))

function validate_config(cfg)
    normalized = static_run_config(cfg)
    require(3normalized["Lx"]*normalized["Ly"] <= 18,
        "example requires a cylinder with at most 18 sites")
    require(normalized["nsweeps"] <= 4, "example permits one to four sweeps")
    require(all(<=(64), normalized["maxdim"]) && normalized["initial_linkdim"] <= 64,
        "example permits dimensions at most 64")
    return cfg
end

function require_threads()
    require(Threads.nthreads() == 1, "run the example with --threads=1")
    BLAS.set_num_threads(1)
end

function write_toml(path, value)
    temporary, io = mktemp(dirname(path))
    try
        TOML.print(io, value; sorted=true)
        close(io)
        Base.Filesystem.rename(temporary, path)
    finally
        isopen(io) && close(io)
        ispath(temporary) && rm(temporary)
    end
end

"""Save a completed small trial, then separately record a direct diagnostic reference."""
function solve(output, config_path)
    require_threads()
    cfg = validate_config(TOML.parsefile(config_path))
    mkpath(dirname(abspath(output)))
    solved = run_static(output, cfg)
    # This example-only reference is not a prerequisite or part of the sealed
    # solve. Its failure cannot relabel a completed optimization/checkpoint.
    path = joinpath(output, "direct-diagnostics.toml")
    write_toml(path, Dict("status" => "running", "role" => "direct_measurement_reference"))
    try
        write_toml(path, static_diagnostics(solved.result))
    catch err
        write_toml(path, Dict("status" => "failed", "role" => "direct_measurement_reference",
            "exception_type" => string(nameof(typeof(err)))))
        rethrow()
    end
    return solved.record
end

"""Diagnose the pinned trial through the shared run loader; preserve all inputs."""
function diagnose(solve_output, output)
    require_threads()
    validate_config(TOML.parsefile(joinpath(solve_output, "config.toml")))
    return diagnose_static_run(solve_output; output)
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
