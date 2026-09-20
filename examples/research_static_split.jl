#!/usr/bin/env julia
# Reuse the existing allowlisted configuration, strict parent load, and physics
# diagnostics. Including the file does not call its guarded main function.
include("research_static.jl")

const SPLIT_BACKEND = "outputs/source-snapshots/n27-refine-before-chirality-20260919"
const SPLIT_SOURCES = ("examples/research_static_split.jl", "examples/run_static_split.py",
    "examples/research_static.jl", "examples/run_static27.py", "examples/run_static18.py")
const SPLIT_LIMITS = Dict("solve"=>900.0,"diagnose"=>180.0)
split_unchanged(hashes) = all(isfile(rootpath(p)) && file_hash(rootpath(p))==h for (p,h) in hashes)

function split_config(path)
    cfg = read_config(path)
    cfg["Lx"]==cfg["Ly"]==3 && cfg["Q"] in (1,5) && cfg["maxdim"]==512 &&
        cfg["batches"]==cfg["batch_sweeps"]==1 && cfg["initialization"]=="resume" ||
        error("only the declared N27 Q1/Q5 one-sweep chi512 cases are supported")
    cfg["initial_completed_sweeps"]==(cfg["Q"]==1 ? 8 : 9) || error("wrong parent sweep count")
    cfg["stationarity"]==Dict("energy_per_site"=>1e-6,"sz_profile"=>1e-4,
        "bond_profile"=>1e-4,"variance_per_site"=>1e-5,"truncation"=>1e-6) ||
        error("the declared precision thresholds must be retained")
    return cfg
end

function split_settings(cfg)
    return (; seed=cfg["seed"],nsweeps=1,maxdim=512,cutoff=0.0,noise=0.0,
        eigsolve_tol=1e-11,eigsolve_krylovdim=40,eigsolve_maxiter=20,measure_variance=false)
end

function split_actual_settings(settings,cfg)
    settings.seed==cfg["seed"] && settings.nsweeps==1 && settings.maxdim==[512] &&
        settings.cutoff==settings.noise==0.0 && settings.eigsolve_tol==1e-11 &&
        settings.eigsolve_krylovdim==40 && settings.eigsolve_maxiter==20 &&
        !settings.measure_variance && settings.initialization=="provided_mps" ||
        error("actual saved solver settings do not match the declared solve")
end

function split_parent!(cfg,lattice,record,inputs)
    saved,previous,hashes = parent_state(cfg,lattice,record)
    merge!(inputs,hashes)
    saved.settings.maxdim==[cfg["Q"]==1 ? 256 : 512] || error("wrong parent chi")
    saved.theta_path==[0.0] && saved.gauge==:seam && all(iszero,saved.hz) || error("wrong parent model/flux")
    previous["energy"]==saved.diagnostics.energy && previous["sz_profile"]==saved.diagnostics.sz ||
        error("parent record and checkpoint observables differ")
    haskey(previous,"bond_energy") || error("parent bond profile is required")
    return saved,previous
end

function split_solve!(output,cfg,lattice,record,row,inputs,persist)
    parent,previous = split_parent!(cfg,lattice,record,inputs)
    row["active_phase"] = "DMRG"
    persist()
    function progress(event)
        selected_bond = event.kind==:bond && event.bond in (1,9,18,nsites(lattice)-1)
        if event.kind==:sweep || event.kind==:phase || selected_bond
            row["last_event"] = Dict(string(k)=>(v isa Symbol ? string(v) : v) for (k,v) in pairs(event))
            persist()
            if selected_bond
                println(cfg["case_id"]," half_sweep=",event.half_sweep," bond=",event.bond,
                    " E=",event.energy," elapsed=",round(elapsed();digits=2))
            else
                println(cfg["case_id"]," phase=",event.kind," elapsed=",round(elapsed();digits=2))
            end
            flush(stdout)
        end
    end
    settings = split_settings(cfg)
    row["dmrg_seconds"] = @elapsed result = run_dmrg(lattice,0.0;Q=cfg["Q"],
        sites=parent.sites,psi0=parent.psi,settings...,progress_callback=progress)
    split_actual_settings(result.settings,cfg)
    length(result.sweep_energies)==length(result.max_truncation_errors)==1 || error("expected exactly one sweep")
    merge!(row,Dict("energy"=>result.energy,"energy_per_site"=>result.energy/nsites(lattice),
        "sz_profile"=>result.sz,"sweep_energies"=>result.sweep_energies,
        "measured_truncation_errors"=>result.max_truncation_errors,
        "maxlinkdim"=>maxlinkdim(result.psi),"settings"=>asdict(result.settings),
        "last_optimizer_energy_error"=>abs(result.energy-only(result.sweep_energies)),
        "within_batch_sweep_change_measured"=>false,"active_phase"=>"checkpoint"))
    persist()
    row["checkpoint_save_seconds"] = @elapsed snapshot = save_checkpoint(output,result;
        baseline=result,theta_path=[0.0],status=:trial)
    row["checkpoint"] = relpath(snapshot,ROOT)
    row["checkpoint_sha256"] = Dict(relpath(joinpath(snapshot,p),ROOT)=>file_hash(joinpath(snapshot,p))
        for p in ("metadata.toml","state.jls","checksums.toml"))
    merge!(inputs,row["checkpoint_sha256"])
    row["checkpoint_metadata_sha256"] = file_hash(joinpath(snapshot,"metadata.toml"))
    row["checkpoint_payload_sha256"] = file_hash(joinpath(snapshot,"state.jls"))
    row["active_phase"] = "strict_checkpoint_reload"
    persist()
    row["checkpoint_reload_seconds"] = @elapsed saved = load_checkpoint(snapshot,lattice;
        Q=cfg["Q"],gauge=:seam,status=:trial,expected_theta=0.0,expected_settings=result.settings)
    row["checkpoint_overlap_error"] = abs(abs(inner(saved.psi,result.psi))-1)
    row["checkpoint_verified"] = row["checkpoint_overlap_error"]<=1e-12 &&
        saved.diagnostics.energy==result.energy && saved.diagnostics.sz==result.sz
    row["checkpoint_verified"] || error("saved trial strict reload mismatch")
    row["status"] = "checkpoint_saved_diagnostics_deferred"
    row["active_phase"] = "completed_solve"
    record["completed_sweeps"] = row["cumulative_sweeps"]
    record["status"] = "completed_solve_diagnostics_deferred"
    # Deliberately no integrity_passed, variance, precision, or stationarity
    # conclusion: this stage establishes only completed solve + strict storage.
    return nothing
end

function split_bind_solve!(solve_output,validation_sha,execution_sha,cfg,lattice,record,inputs)
    validation_path,execution_path = joinpath(solve_output,"validation.toml"),joinpath(solve_output,"execution.toml")
    inputs[relpath(validation_path,ROOT)] = validation_sha
    inputs[relpath(execution_path,ROOT)] = execution_sha
    split_unchanged(inputs) || error("solve record or execution pin mismatch")
    old,execution = TOML.parsefile(validation_path),TOML.parsefile(execution_path)
    old["stage"]=="solve" && old["status"]=="completed_solve_diagnostics_deferred" &&
        execution["status"]=="exited" && execution["worker_exit_confirmed"]===true &&
        execution["worker_exit_code"]==0 || error("solve has not completed successfully")
    all(old[k]===true for k in ("sources_unchanged","backend_unchanged","inputs_unchanged","parent_unchanged")) ||
        error("solve final provenance guards did not pass")
    old["config"]==cfg && old["config_sha256"]==record["config_sha256"] &&
        old["configuration"]==record["configuration"] && old["runtime"]==record["runtime"] &&
        old["code"]["source_sha256"]==record["code"]["source_sha256"] &&
        old["solver"]==record["solver"] || error("solve configuration/runtime/backend mismatch")
    command = execution["command"]
    length(command)==8 && basename(command[1])=="julia" && command[2:4]==
        ["--project=.","--startup-file=no","--threads=1"] &&
        rootpath(command[5])==joinpath(ROOT,"examples/research_static_split.jl") && command[6]=="solve" &&
        rootpath(command[7])==solve_output && rootpath(command[8])==rootpath(old["config_path"]) ||
        error("solve execution command does not bind this case")
    execution["wall_limit_seconds"]==old["wall_limit_seconds"] &&
        0<execution["wall_limit_seconds"]<=SPLIT_LIMITS["solve"] || error("invalid solve wall budget")
    for (p,h) in old["analysis_source_sha256"]
        inputs[relpath(joinpath(solve_output,"analysis-sources",p),ROOT)] = h
    end
    inputs[relpath(joinpath(solve_output,"analysis-sources/config.toml"),ROOT)] = old["config_sha256"]
    merge!(inputs,old["input_sha256"])
    original = only(old["batches"])
    original["status"]=="checkpoint_saved_diagnostics_deferred" && original["checkpoint_verified"]===true &&
        original["cumulative_sweeps"]==cfg["initial_completed_sweeps"]+1==old["completed_sweeps"] &&
        original["within_batch_sweep_change_measured"]===false || error("wrong solve row")
    all(!haskey(original,k) for k in ("integrity_passed","variance","last_sweep_energy_change")) ||
        error("unexpected solve diagnostic claim")
    snapshot = rootpath(original["checkpoint"])
    dirname(dirname(snapshot))==solve_output || error("checkpoint is not owned by solve output")
    expected_files = Set(relpath(joinpath(snapshot,p),ROOT) for p in ("metadata.toml","state.jls","checksums.toml"))
    Set(keys(original["checkpoint_sha256"]))==expected_files || error("incomplete checkpoint pins")
    merge!(inputs,original["checkpoint_sha256"])
    split_unchanged(inputs) || error("solve archive/checkpoint/parent bytes changed")
    meta = TOML.parsefile(joinpath(snapshot,"metadata.toml"))
    meta["settings"]==original["settings"] && meta["configuration"]==record["configuration"] &&
        meta["provenance"]["source_sha256"]==record["code"]["source_sha256"] &&
        merge(meta["runtime"],Dict("julia_threads"=>1,"blas_threads"=>1))==record["runtime"] ||
        error("checkpoint metadata does not bind solve record")
    settings = KagomeDMRG._checkpoint_settings(meta["settings"])
    split_actual_settings(settings,cfg)
    state = meta["state"]
    state["energy"]==original["energy"] && state["sz"]==original["sz_profile"] &&
        state["sweep_energies"]==original["sweep_energies"] &&
        state["max_truncation_errors"]==original["measured_truncation_errors"] || error("solve diagnostics mismatch")
    record["solve_record"] = relpath(validation_path,ROOT)
    record["solve_record_sha256"] = validation_sha
    record["solve_execution"] = relpath(execution_path,ROOT)
    record["solve_execution_sha256"] = execution_sha
    return original,snapshot,settings
end

function split_diagnose!(solve_output,validation_sha,execution_sha,cfg,lattice,record,row,inputs,persist)
    original,snapshot,settings = split_bind_solve!(solve_output,validation_sha,execution_sha,cfg,lattice,record,inputs)
    parent,previous = split_parent!(cfg,lattice,record,inputs)
    merge!(row,deepcopy(original))
    row["status"],row["active_phase"] = "running","strict_checkpoint_load"
    persist()
    row["diagnostic_checkpoint_load_seconds"] = @elapsed saved = load_checkpoint(snapshot,lattice;
        Q=cfg["Q"],gauge=:seam,status=:trial,expected_theta=0.0,expected_settings=settings)
    saved.theta_path==[0.0] && saved.diagnostics.energy==row["energy"] &&
        saved.diagnostics.sz==row["sz_profile"] || error("strictly loaded state differs from solve record")
    row["active_phase"] = "remeasure_energy_density_norm_charge"
    persist()
    H = twisted_exchange_mpo(saved.sites,lattice,0.0;gauge=:seam,hz=saved.hz)
    energy,sz = inner(saved.psi',H,saved.psi),sz_profile(saved.psi)
    row["remeasured_energy_real"],row["remeasured_energy_imag"] = real(energy),imag(energy)
    row["remeasured_energy_error"] = abs(energy-row["energy"])
    row["remeasured_sz_error"] = maximum(abs.(sz-row["sz_profile"]))
    row["remeasured_norm_error"] = abs(norm(saved.psi)-1)
    row["remeasured_charge_error"] = abs(sum(sz)-cfg["Q"]/2)
    row["saved_state_remeasurement_passed"] = isfinite(energy) && all(isfinite,sz) &&
        row["remeasured_energy_error"]<=1e-10*max(1,abs(row["energy"])) &&
        row["remeasured_sz_error"]<=1e-12 && row["remeasured_norm_error"]<=1e-12 &&
        row["remeasured_charge_error"]<=1e-12
    row["saved_state_remeasurement_passed"] || error("saved state remeasurement failed")
    row["parent_overlap_abs"] = abs(inner(parent.psi,saved.psi))
    row["parent_overlap_measured"] = true
    row["column_sz"] = [sum(sz[3lattice.Ly*x+1:3lattice.Ly*(x+1)]) for x in 0:lattice.Lx-1]
    # Keep the original energy/profile convention; fresh checks above establish
    # agreement, while no diagnostic rewrites the checkpoint or its settings.
    result = (; psi=saved.psi,H,Q=cfg["Q"],energy=row["energy"],sz=row["sz_profile"],
        max_truncation_errors=saved.diagnostics.max_truncation_errors)
    measure_state(result,lattice,record,row,persist)
    row["energy_change_from_previous"] = row["energy"]-previous["energy"]
    row["max_sz_change"] = maximum(abs.(row["sz_profile"]-previous["sz_profile"]))
    row["max_bond_change"] = maximum(abs.(row["bond_energy"]-previous["bond_energy"]))
    values = Dict("energy_per_site"=>max(abs(row["energy_change_from_previous"]),row["last_optimizer_energy_error"])/nsites(lattice),
        "sz_profile"=>row["max_sz_change"],"bond_profile"=>row["max_bond_change"],
        "variance_per_site"=>abs(row["variance_per_site"]),"truncation"=>only(row["measured_truncation_errors"]))
    row["precision_values"] = values
    row["precision_passed"] = Dict(k=>isfinite(v) && v<=cfg["stationarity"][k] for (k,v) in values)
    row["all_precision_conditions_passed"] = all(Base.values(row["precision_passed"]))
    row["same_chi_comparison"] = parent.settings.maxdim==settings.maxdim
    row["comparison_kind"] = row["same_chi_comparison"] ? "same_chi" : "chi_change"
    row["stationarity_evaluated"] = row["same_chi_comparison"]
    row["stationarity_passed"] = row["same_chi_comparison"] && row["all_precision_conditions_passed"]
    row["integrity_passed"] || error("observable consistency failed")
    record["completed_sweeps"] = row["cumulative_sweeps"]
    record["status"] = "completed_saved_state_diagnostics_accuracy_separate"
    return nothing
end

function split_main(stage,output,config_path,extra)
    haskey(SPLIT_LIMITS,stage) || error("expected solve or diagnose")
    isdir(output) && readdir(output)==["execution.toml"] || error("use fresh launcher output")
    execution = TOML.parsefile(joinpath(output,"execution.toml"))
    execution["status"]=="running" && 0<execution["wall_limit_seconds"]<=SPLIT_LIMITS[stage] || error("invalid supervisor")
    Threads.nthreads()==BLAS.get_num_threads()==1 || error("one Julia/BLAS thread required")
    realpath(pkgdir(KagomeDMRG))==realpath(rootpath(SPLIT_BACKEND)) || error("original archived backend required")
    cfg = split_config(config_path)
    sources = Dict(p=>file_hash(rootpath(p)) for p in SPLIT_SOURCES)
    inputs = Dict(relpath(config_path,ROOT)=>file_hash(config_path))
    lattice = kagome_cylinder(cfg["Lx"],cfg["Ly"])
    row = Dict{String,Any}("batch"=>1,"status"=>"running","active_phase"=>"verify_parent",
        "cumulative_sweeps"=>cfg["initial_completed_sweeps"]+1)
    record = Dict{String,Any}("schema_version"=>1,"stage"=>stage,"status"=>"running",
        "case_id"=>cfg["case_id"],"recorded_at_utc"=>string(now(UTC)),"N"=>nsites(lattice),"Q"=>cfg["Q"],"theta"=>0.0,
        "configuration"=>KagomeDMRG._checkpoint_configuration(lattice,:seam,nothing;Q=cfg["Q"]),
        "config"=>cfg,"config_path"=>relpath(config_path,ROOT),"config_sha256"=>file_hash(config_path),
        "analysis_source_sha256"=>sources,"input_sha256"=>inputs,"code"=>KagomeDMRG._checkpoint_provenance(),
        "backend_project"=>SPLIT_BACKEND,"backend_snapshot_origin"=>cfg["backend_snapshot_origin"],
        "git_revision_semantics"=>"execution_directory_git_lookup_not_snapshot_origin",
        "runtime"=>merge(KagomeDMRG._checkpoint_runtime(),Dict("julia_threads"=>1,"blas_threads"=>1)),
        "solver"=>asdict(split_settings(cfg)),"integrity_limits"=>INTEGRITY_LIMITS,
        "wall_limit_seconds"=>execution["wall_limit_seconds"],"new_DMRG"=>stage=="solve","new_ED"=>false,
        "diagnostics_deferred"=>stage=="solve","solve_records_modified"=>false,
        "original_solve_global_integrity_retroactively_claimed"=>false,
        "phase_identification"=>"not_attempted","bulk_plateau"=>"not_established",
        "variance_interpretation"=>"H_dispersion_not_ground_energy_error_bound",
        "stationarity_interpretation"=>"fixed_chi_observation_requires_chi_seed_size_checks",
        "batches"=>[row])
    persist() = begin
        record["worker_elapsed_seconds"] = elapsed()
        KagomeDMRG._write_flux_record(joinpath(output,"validation.toml"),record)
    end
    persist()
    try
        for p in SPLIT_SOURCES
            destination = joinpath(output,"analysis-sources",p)
            mkpath(dirname(destination))
            cp(rootpath(p),destination)
        end
        cp(config_path,joinpath(output,"analysis-sources/config.toml"))
        if stage=="solve"
            isempty(extra) || error("solve takes no diagnostic inputs")
            split_solve!(output,cfg,lattice,record,row,inputs,persist)
        else
            length(extra)==3 || error("diagnose requires SOLVE_OUTPUT VALIDATION_SHA EXECUTION_SHA")
            split_diagnose!(abspath(extra[1]),extra[2],extra[3],cfg,lattice,record,row,inputs,persist)
        end
    catch exception
        record["status"] = exception isa InterruptException ? "interrupted" : "failed"
        record["exception_type"] = string(nameof(typeof(exception)))
        rethrow()
    finally
        try
            record["sources_unchanged"] = split_unchanged(sources) && all(
                file_hash(joinpath(output,"analysis-sources",p))==h for (p,h) in sources) &&
                file_hash(joinpath(output,"analysis-sources/config.toml"))==record["config_sha256"]
            record["inputs_unchanged"] = split_unchanged(inputs)
            record["parent_unchanged"] = haskey(record,"parent") && split_unchanged(record["parent"]["checkpoint_sha256"])
            record["backend_unchanged"] = KagomeDMRG._checkpoint_provenance()["source_sha256"]==record["code"]["source_sha256"]
            all(record[k] for k in ("sources_unchanged","inputs_unchanged","parent_unchanged","backend_unchanged")) ||
                (record["status"]="source_or_input_changed")
        catch exception
            record["status"] = "source_or_input_check_failed"
            record["final_check_exception_type"] = string(nameof(typeof(exception)))
        end
        persist()
    end
    return record
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)>=3 || error("usage: research_static_split.jl STAGE OUTPUT CONFIG [SOLVE_OUTPUT VALIDATION_SHA EXECUTION_SHA]")
    record = split_main(ARGS[1],abspath(ARGS[2]),abspath(ARGS[3]),ARGS[4:end])
    println(record["status"])
    record["status"] in ("completed_solve_diagnostics_deferred","completed_saved_state_diagnostics_accuracy_separate") || exit(1)
end
