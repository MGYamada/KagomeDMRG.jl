const _STATIC_RUN_FORMAT = "KagomeDMRG.static_run"
const _STATIC_RUN_REQUIRED_CONFIG = ("schema_version", "Lx", "Ly", "Q", "seed", "nsweeps", "maxdim")
const _STATIC_RUN_DEFAULTS = (
    Jxy=1.0, Jz=1.0, gauge="seam", initial_linkdim=4, cutoff=1e-12,
    noise=0.0, eigsolve_tol=1e-12, eigsolve_krylovdim=20, eigsolve_maxiter=10)
const _STATIC_RUN_FILES = ("metadata.toml", "state.jls", "checksums.toml")

_static_run_require(condition, message) = condition || throw(ArgumentError(message))

# TOML integer/float/bool values can compare equal in Julia. A saved contract
# must retain their roles as well as their values; dictionary/vector container
# element types need not survive TOML parsing.
function _static_run_matches(actual, expected)
    if expected isa AbstractDict
        return actual isa AbstractDict && Set(keys(actual)) == Set(keys(expected)) &&
            all(k -> _static_run_matches(actual[k], expected[k]), keys(expected))
    elseif expected isa AbstractVector
        return actual isa AbstractVector && length(actual) == length(expected) &&
            all(_static_run_matches(a,b) for (a,b) in zip(actual,expected))
    end
    return typeof(actual) === typeof(expected) && isequal(actual, expected)
end

function _static_run_int(value, name; minimum=typemin(Int))
    _static_run_require(value isa Integer && !(value isa Bool) &&
        minimum <= value <= typemax(Int), "$name must be an integer in range")
    return Int(value)
end

function _static_run_real(value, name; minimum=-Inf, positive=false)
    _static_run_require(value isa Real && !(value isa Bool), "$name must be a real number")
    number = Float64(value)
    _static_run_require(isfinite(number) && number >= minimum && (!positive || number > 0),
        "$name must be finite and in range")
    return number
end

"""
    static_run_config(config::AbstractDict)

Validate and copy a string-keyed, schema-1 static-run configuration. Required
keys are `schema_version`, `Lx`, `Ly`, `Q`, `seed`, `nsweeps`, and `maxdim`.
Optional keys are `Jxy`, `Jz`, `hz`, `gauge`, `initial_linkdim`, `cutoff`, `noise`,
`eigsolve_tol`, `eigsolve_krylovdim`, and `eigsolve_maxiter`. Return a detached
TOML-compatible dictionary with all defaults materialized. Reject unknown keys,
Bool-valued numbers, nonfinite values, invalid charge, and schedules longer than
`nsweeps`. `maxdim` is a vector; its last entry repeats over remaining sweeps.

This workflow prepares a zero-flux nearest-neighbor trial from a random complex
fixed-charge MPS and defers variance to diagnostics. It imposes no research
size/sweep/dimension budget; the small example imposes its own limits.
"""
function static_run_config(config::AbstractDict)
    allowed = Set(vcat(collect(_STATIC_RUN_REQUIRED_CONFIG),
        string.(collect(keys(_STATIC_RUN_DEFAULTS))), ["hz"]))
    _static_run_require(all(k -> k isa AbstractString, keys(config)) &&
        issubset(Set(keys(config)), allowed) &&
        all(k -> haskey(config, k), _STATIC_RUN_REQUIRED_CONFIG),
        "configuration has missing or unknown keys")
    cfg = merge(Dict{String,Any}(string(k) => v for (k,v) in pairs(_STATIC_RUN_DEFAULTS)),
                Dict{String,Any}(config))
    for name in ("schema_version", "Lx", "Ly", "seed", "nsweeps", "initial_linkdim",
                 "eigsolve_krylovdim", "eigsolve_maxiter")
        minimum = name == "seed" ? 0 : name == "Ly" ? 3 : name == "eigsolve_krylovdim" ? 2 : 1
        cfg[name] = _static_run_int(cfg[name], name; minimum)
    end
    _static_run_require(cfg["schema_version"] == 1, "unsupported static configuration schema")
    N = try
        Base.checked_mul(Base.checked_mul(3, cfg["Lx"]), cfg["Ly"])
    catch err
        err isa OverflowError || rethrow()
        throw(ArgumentError("static configuration site count overflows Int"))
    end
    cfg["Q"] = _static_run_int(cfg["Q"], "Q")
    target_sector(N; Q=cfg["Q"])
    dims = cfg["maxdim"]
    _static_run_require(dims isa AbstractVector && 1 <= length(dims) <= cfg["nsweeps"],
        "maxdim must be a nonempty vector no longer than nsweeps")
    cfg["maxdim"] = [_static_run_int(d, "maxdim"; minimum=1) for d in dims]
    for name in ("Jxy", "Jz", "cutoff", "noise", "eigsolve_tol")
        cfg[name] = _static_run_real(cfg[name], name;
            minimum=name in ("Jxy", "Jz") ? -Inf : 0.0, positive=name == "eigsolve_tol")
    end
    _static_run_require(cfg["gauge"] isa AbstractString && cfg["gauge"] in ("seam", "uniform"),
        "gauge must be seam or uniform")
    cfg["gauge"] = String(cfg["gauge"])
    fields = get(cfg, "hz", nothing)
    if fields === nothing && !haskey(config, "hz")
        fields = zeros(N)
    end
    _static_run_require(fields isa AbstractVector && length(fields) == N,
        "hz must be a vector with one field per site")
    cfg["hz"] = [_static_run_real(h, "hz") for h in fields]
    return cfg
end

_static_run_lattice(cfg) = kagome_cylinder(cfg["Lx"], cfg["Ly"]; Jxy=cfg["Jxy"], Jz=cfg["Jz"])
_static_run_settings(cfg) = (; seed=cfg["seed"], initial_linkdim=cfg["initial_linkdim"],
    nsweeps=cfg["nsweeps"], maxdim=copy(cfg["maxdim"]), cutoff=cfg["cutoff"], noise=cfg["noise"],
    eigsolve_tol=cfg["eigsolve_tol"], eigsolve_krylovdim=cfg["eigsolve_krylovdim"],
    eigsolve_maxiter=cfg["eigsolve_maxiter"], measure_variance=false)

# Resolve the existing parent before testing containment, including symlink aliases.
function _static_run_destination(output; outside=nothing)
    absolute = abspath(output)
    destination = joinpath(realpath(dirname(absolute)), basename(absolute))
    _static_run_require(!ispath(destination) && !islink(destination), "output must be a new directory")
    outside === nothing || _static_run_require(
        first(splitpath(relpath(destination, outside))) == "..",
        "diagnostic output must be outside the solve directory")
    return destination
end

function _static_run_write_seal(path, hash)
    temporary, io = mktemp(dirname(path))
    try
        write(io, hash * "\n")
        close(io)
        Base.Filesystem.rename(temporary, path)
    finally
        isopen(io) && close(io)
        ispath(temporary) && rm(temporary)
    end
end

"""
    run_static(output, config; progress_callback=nothing)

Optimize the validated zero-flux static configuration and save a trial checkpoint
before any optional diagnostic measurement. `output` must be new and its parent
must exist. Atomically update `solve.toml` through running/failed/completed_solve,
then publish `solve.sha256` only on success. Return `(record, result, snapshot)`;
these in-memory values may be changed by the caller without rewriting saved data.

`diagnostic_status` remains `not_run`. Completion is not convergence or branch
acceptance. A callback has the `run_dmrg` event contract; exceptions propagate,
with only their type and the active phase saved. Interrupted work may retain a
running record or an unsealed completed record, neither of which can be loaded
as a completed run. The library records but does not set Julia/BLAS threads.
"""
function run_static(output::AbstractString, config::AbstractDict; progress_callback=nothing)
    cfg = static_run_config(config)
    identity = _execution_identity()
    lattice = _static_run_lattice(cfg)
    configuration = _checkpoint_configuration(lattice, Symbol(cfg["gauge"]), cfg["hz"]; Q=cfg["Q"])
    destination = _static_run_destination(output)
    record = Dict{String,Any}("schema_version" => 1, "format" => _STATIC_RUN_FORMAT,
        "stage" => "solve", "status" => "running", "active_phase" => "solve",
        "solver_phase" => "not_started", "diagnostic_status" => "not_run",
        "config" => deepcopy(cfg), "configuration" => configuration,
        "theta" => 0.0, "checkpoint_status" => "trial",
        "settings" => _static_dict(merge(_static_run_settings(cfg), (; initialization="random_fixed_charge"))),
        "runtime" => Dict(identity.runtime), "provenance" => _checkpoint_provenance(identity),
        "threads" => Dict("julia" => Threads.nthreads(), "blas" => BLAS.get_num_threads()),
        "timing_seconds" => Dict{String,Float64}())
    mkdir(destination)
    record_path = joinpath(destination, "solve.toml")
    persist() = _write_flux_record(record_path, record)
    started = time_ns()
    try
        _write_flux_record(joinpath(destination, "config.toml"), cfg)
        persist()
        function progress(event)
            if event.kind == :phase
                record["solver_phase"] = string(event.phase)
                persist()
            end
            progress_callback === nothing || progress_callback(event)
            return nothing
        end
        record["timing_seconds"]["solve"] = @elapsed result = run_dmrg(lattice, 0.0;
            Q=cfg["Q"], gauge=Symbol(cfg["gauge"]), hz=cfg["hz"],
            _static_run_settings(cfg)..., progress_callback=progress)
        record["active_phase"] = "checkpoint_save"
        persist()
        record["timing_seconds"]["checkpoint_save"] = @elapsed snapshot = save_checkpoint(
            destination, result; baseline=result, theta_path=[0.0], status=:trial)
        record["checkpoint"] = relpath(snapshot, destination)
        names = vcat(["config.toml"], [joinpath(record["checkpoint"], n) for n in _STATIC_RUN_FILES])
        record["file_sha256"] = Dict(n => _file_sha256(joinpath(destination,n)) for n in names)
        _static_run_require(identity == _execution_identity(), "implementation changed during static run")
        record["timing_seconds"]["total"] = Float64(time_ns()-started)*1e-9
        record["status"], record["active_phase"] = "completed_solve", "completed"
        persist()
        _static_run_write_seal(joinpath(destination, "solve.sha256"), _file_sha256(record_path))
        return (; record=deepcopy(record), result, snapshot)
    catch err
        record["status"] = "failed"
        record["exception_type"] = string(nameof(typeof(err)))
        record["timing_seconds"]["total"] = Float64(time_ns()-started)*1e-9
        persist()
        rethrow()
    end
end

function _static_run_unchanged(pins)
    _static_run_require(all(isfile(path) && _file_sha256(path) == hash for (path,hash) in pins),
        "a static run input or original solve record changed")
end

function _read_static_run(directory)
    identity = _execution_identity()
    root = realpath(directory)
    record_path, seal = joinpath(root,"solve.toml"), joinpath(root,"solve.sha256")
    record_hash = _file_sha256(record_path)
    _static_run_require(strip(read(seal, String)) == record_hash, "solve record checksum mismatch")
    record = TOML.parsefile(record_path)
    expected_keys = Set(["schema_version", "format", "stage", "status", "active_phase",
        "solver_phase", "diagnostic_status", "config", "configuration", "theta", "checkpoint_status",
        "settings", "runtime", "provenance", "threads", "timing_seconds", "checkpoint", "file_sha256"])
    _static_run_require(Set(keys(record)) == expected_keys &&
        record["schema_version"] isa Integer && !(record["schema_version"] isa Bool) &&
        record["schema_version"] == 1 && record["format"] == _STATIC_RUN_FORMAT,
        "unsupported static run schema or format")
    _static_run_require(record["stage"] == "solve" && record["status"] == "completed_solve" &&
        record["active_phase"] == "completed" && record["solver_phase"] == "sz" &&
        record["diagnostic_status"] == "not_run" && record["theta"] === 0.0 &&
        record["checkpoint_status"] == "trial", "expected a completed zero-flux trial solve")
    cfg = static_run_config(record["config"])
    _static_run_require(_static_run_matches(record["config"], cfg) &&
        _static_run_matches(TOML.parsefile(joinpath(root,"config.toml")), cfg),
        "solve configuration mismatch")
    lattice = _static_run_lattice(cfg)
    configuration = _checkpoint_configuration(lattice, Symbol(cfg["gauge"]), cfg["hz"]; Q=cfg["Q"])
    settings = _static_dict(merge(_static_run_settings(cfg), (; initialization="random_fixed_charge")))
    _static_run_require(_static_run_matches(record["configuration"], configuration) &&
        _static_run_matches(record["settings"], settings),
        "solve model or settings do not match configuration")
    provenance = record["provenance"]
    _static_run_require(_static_run_matches(record["runtime"], Dict(identity.runtime)) &&
        provenance["source_sha256"] == Dict(identity.source_sha256) &&
        provenance["environment_sha256"] == Dict(identity.environment_sha256) &&
        provenance["manifest"] == identity.manifest, "static run source or environment mismatch")
    _static_run_require(Set(keys(record["threads"])) == Set(["julia", "blas"]) &&
        all(v -> v isa Integer && !(v isa Bool) && v > 0, values(record["threads"])) &&
        Set(keys(record["timing_seconds"])) == Set(["solve", "checkpoint_save", "total"]) &&
        all(v -> v isa Real && !(v isa Bool) && isfinite(v) && v >= 0, values(record["timing_seconds"])),
        "invalid static run timing or thread metadata")
    parts = splitpath(record["checkpoint"])
    _static_run_require(length(parts) == 2 && parts[1] == "trial" && startswith(parts[2], "checkpoint-"),
        "checkpoint must belong to the solve output")
    snapshot = joinpath(root, record["checkpoint"])
    names = Set(vcat(["config.toml"], [joinpath(record["checkpoint"], n) for n in _STATIC_RUN_FILES]))
    _static_run_require(Set(keys(record["file_sha256"])) == names, "incomplete solve artifact pins")
    _static_run_require(all(h -> h isa String && occursin(r"^[0-9a-f]{64}$", h),
        values(record["file_sha256"])), "invalid solve artifact hash")
    # Paths are canonical local files, never aliases to another run's artifacts.
    for name in union(names, Set(["solve.toml", "solve.sha256"]))
        path = joinpath(root, name)
        _static_run_require(realpath(path) == path, "static run artifacts must not traverse symlinks")
    end
    pins = Dict(joinpath(root,n) => hash for (n,hash) in record["file_sha256"])
    pins[record_path], pins[seal] = record_hash, _file_sha256(seal)
    _static_run_unchanged(pins)
    return (; root, config=cfg, record, snapshot, lattice, pins, record_sha256=record_hash)
end

function _checked_static_run(directory)
    try
        return _read_static_run(directory)
    catch err
        err isa InterruptException && rethrow()
        err isa ArgumentError && rethrow()
        throw(ArgumentError("invalid static run ($(nameof(typeof(err))))"))
    end
end

"""
    load_static_run(directory)

Validate a sealed schema-1 run and strictly load its trial checkpoint. Return
`(config, record, snapshot, checkpoint, record_sha256)`. Check the full model,
settings, artifact allowlist, checksums, source and active environment before
deserializing. Artifact paths cannot traverse symlinks. Inputs remain unchanged.
Only trusted local files from the same source/runtime/environment are supported;
older ad hoc example records are not migrated or silently reinterpreted.
"""
function load_static_run(directory::AbstractString)
    run = _checked_static_run(directory)
    cfg = run.config
    try
        checkpoint = load_checkpoint(run.snapshot, run.lattice; Q=cfg["Q"],
            gauge=Symbol(cfg["gauge"]), hz=cfg["hz"], expected_theta=0.0,
            expected_settings=run.record["settings"], status=:trial)
        return (; config=cfg, record=run.record, snapshot=run.snapshot, checkpoint,
                record_sha256=run.record_sha256)
    finally
        _static_run_unchanged(run.pins)
    end
end

"""
    diagnose_static_run(directory; output, cuts=..., chirality=false, reference=nothing)

Validate a completed static run, then diagnose its trial checkpoint in a fresh
directory outside the entire solve directory. Options and returned report follow
`diagnose_checkpoint`. Run-envelope and output-path preflight failures create
no output. `driver.toml` links
the diagnostic to the sealed solve record and reports input preservation, even
on measurement failure. Never update the original solve record or checkpoint.
"""
function diagnose_static_run(directory::AbstractString; output::AbstractString,
        cuts=nothing, chirality::Bool=false, reference=nothing)
    run = _checked_static_run(directory)
    destination = _static_run_destination(output; outside=run.root)
    cfg = run.config
    selected_cuts = cuts === nothing ? (1:run.lattice.Lx-1) : cuts
    try
        return diagnose_checkpoint(run.snapshot, run.lattice; output=destination,
            Q=cfg["Q"], gauge=Symbol(cfg["gauge"]), hz=cfg["hz"], expected_theta=0.0,
            expected_settings=run.record["settings"], status=:trial,
            cuts=selected_cuts, chirality, reference)
    finally
        unchanged = try
            _static_run_unchanged(run.pins)
            true
        catch
            false
        end
        if isdir(destination)
            _write_flux_record(joinpath(destination,"driver.toml"), Dict(
                "schema_version" => 1, "format" => "KagomeDMRG.static_run_diagnostic",
                "solve_record_sha256" => run.record_sha256,
                "solve_artifacts_unchanged" => unchanged,
                "original_solve_records_modified" => !unchanged))
        end
        _static_run_require(unchanged, "static run inputs changed during diagnostics")
    end
end
