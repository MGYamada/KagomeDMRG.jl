#!/usr/bin/env julia
# Bounded J1=1,J2=J3=0.5,Q=0 measurement-control work. The launcher owns
# execution.toml and the inclusive timeout; this worker owns validation.toml.
const CSL_WORKER_START = time_ns()
using KagomeDMRG, ITensors, ITensorMPS, LinearAlgebra, SHA, TOML, Dates
BLAS.set_num_threads(1)
ITensors.disable_threaded_blocksparse()

const CSL_ROOT = normpath(joinpath(@__DIR__, ".."))
const CSL_COUPLINGS = (; J1=1.0, J2=0.5, J3=0.5)
const CSL_CUTS = [1, 2]
const CSL_SCHMIDT_BONDS = [12, 24]
const CSL_CENTER_SITES = collect(13:24)
const CSL_EXTRA_SOURCES = ("examples/research_csl.jl", "examples/run_research.py",
    "examples/run_static18.py", "examples/run_static27.py")
const CSL_SOLVER_DEFAULTS = Dict{String,Any}("seed"=>11, "maxdim"=>128,
    "initial_linkdim"=>4, "eigsolve_tol"=>1e-11,
    "eigsolve_krylovdim"=>40, "eigsolve_maxiter"=>20)
const CSL_INTEGRITY_LIMITS = Dict("state_norm_error"=>1e-10,
    "total_sz_error"=>1e-9, "scaled_bond_energy_sum_error"=>1e-9,
    "schmidt_density_error"=>1e-9, "schmidt_probability_sum_error"=>1e-10,
    "checkpoint_overlap_error"=>1e-10)

csl_elapsed() = Float64(time_ns()-CSL_WORKER_START)/1e9
csl_dict(value) = Dict{String,Any}(string(k)=>v for (k,v) in pairs(value))
csl_sha(path) = bytes2hex(sha256(read(path)))
csl_extra_hashes() = Dict(path=>csl_sha(joinpath(CSL_ROOT,path)) for path in CSL_EXTRA_SOURCES)
csl_path(path) = isabspath(path) ? normpath(path) : normpath(joinpath(CSL_ROOT,path))
csl_record_path(path) = relpath(abspath(path), CSL_ROOT)

function csl_keys(table, allowed, label)
    table isa AbstractDict || throw(ArgumentError("$label must be a table"))
    isempty(setdiff(Set(keys(table)), Set(allowed))) ||
        throw(ArgumentError("unknown key in $label"))
    return table
end

function csl_integer(value, low, high, label)
    value isa Integer && !(value isa Bool) && low <= value <= high ||
        throw(ArgumentError("invalid $label"))
    return Int(value)
end

function csl_finite(value, label; positive=false)
    value isa Real && !(value isa Bool) && isfinite(value) &&
        (positive ? value > 0 : value >= 0) || throw(ArgumentError("invalid $label"))
    return Float64(value)
end

"""Read only the worker's explicit configuration allowlist; paths are ROOT-relative."""
function csl_read_config(path)
    input = TOML.parsefile(path)
    csl_keys(input, ("schema_version", "mode", "start_checkpoint", "start_status",
        "solver", "preparation", "flux", "policy"), "configuration")
    csl_integer(get(input,"schema_version",0),1,1,"schema_version")
    mode = get(input,"mode","")
    mode in ("prepare","flux") || throw(ArgumentError("mode must be prepare or flux"))
    start = get(input,"start_checkpoint","")
    start isa String || throw(ArgumentError("start_checkpoint must be a string"))
    start_status = get(input,"start_status","trial")
    start_status in ("trial","accepted") || throw(ArgumentError("invalid start_status"))
    isempty(start) && start_status != "trial" &&
        throw(ArgumentError("start_status requires a checkpoint"))
    solver_input = get(input,"solver",Dict{String,Any}())
    csl_keys(solver_input,keys(CSL_SOLVER_DEFAULTS),"solver")
    solver = merge(CSL_SOLVER_DEFAULTS,solver_input)
    solver["seed"] = csl_integer(solver["seed"],0,typemax(Int),"seed")
    solver["maxdim"] = csl_integer(solver["maxdim"],2,512,"maxdim")
    solver["initial_linkdim"] = csl_integer(solver["initial_linkdim"],1,128,"initial_linkdim")
    solver["eigsolve_tol"] = csl_finite(solver["eigsolve_tol"],"eigsolve_tol";positive=true)
    solver["eigsolve_krylovdim"] = csl_integer(solver["eigsolve_krylovdim"],2,128,"eigsolve_krylovdim")
    solver["eigsolve_maxiter"] = csl_integer(solver["eigsolve_maxiter"],1,100,"eigsolve_maxiter")
    cfg = Dict{String,Any}("schema_version"=>1,"mode"=>mode,"solver"=>solver,
        "start_checkpoint"=>isempty(start) ? "" : csl_record_path(csl_path(start)),
        "start_status"=>start_status)
    if mode == "prepare"
        !haskey(input,"flux") && !haskey(input,"policy") ||
            throw(ArgumentError("flux configuration is not used by prepare"))
        prep = get(input,"preparation",Dict{String,Any}())
        csl_keys(prep,("batches","batch_sweeps"),"preparation")
        cfg["preparation"] = Dict("batches"=>csl_integer(get(prep,"batches",2),1,2,"batches"),
            "batch_sweeps"=>csl_integer(get(prep,"batch_sweeps",2),2,2,"batch_sweeps"))
    else
        isempty(start) && throw(ArgumentError("flux requires an explicit checkpoint"))
        !haskey(input,"preparation") || throw(ArgumentError("preparation table is not used by flux"))
        flux_input = get(input,"flux",Dict{String,Any}())
        csl_keys(flux_input,("targets","initial_step","min_step","max_trials","sweeps",
            "max_bond_change","max_chirality_change"),"flux")
        required = ("targets","initial_step","min_step","max_trials",
            "max_bond_change","max_chirality_change")
        all(k->haskey(flux_input,k),required) || throw(ArgumentError("missing explicit flux limits"))
        targets = flux_input["targets"]
        targets isa AbstractVector && 1 <= length(targets) <= 32 &&
            all(x->x isa Real && !(x isa Bool) && isfinite(x),targets) ||
            throw(ArgumentError("targets must be a finite nonempty vector"))
        f = Dict{String,Any}("targets"=>Float64.(targets),
            "sweeps"=>csl_integer(get(flux_input,"sweeps",2),2,8,"flux sweeps"),
            "max_trials"=>csl_integer(flux_input["max_trials"],1,64,"max_trials"))
        for key in ("initial_step","min_step","max_bond_change","max_chirality_change")
            f[key] = csl_finite(flux_input[key],key;positive=key in ("initial_step","min_step"))
        end
        f["min_step"] <= f["initial_step"] || throw(ArgumentError("min_step exceeds initial_step"))
        policy_input = get(input,"policy",Dict{String,Any}())
        pkeys = string.(fieldnames(FluxPolicy))
        csl_keys(policy_input,pkeys,"policy")
        all(k->haskey(policy_input,k),pkeys) || throw(ArgumentError("all FluxPolicy limits are explicit"))
        p = Dict(k=>csl_finite(policy_input[k],k;positive=k=="consistency_tol") for k in pkeys)
        FluxPolicy(; (Symbol(k)=>v for (k,v) in p)...)
        cfg["flux"], cfg["policy"] = f, p
    end
    return cfg
end

function csl_snapshot_record(path)
    # load_checkpoint has already verified these hashes before this record is used.
    checksums = TOML.parsefile(joinpath(path,"checksums.toml"))
    return Dict("path"=>csl_record_path(path),
        "checksums_file_sha256"=>csl_sha(joinpath(path,"checksums.toml")),
        "metadata_sha256"=>checksums["files"]["metadata.toml"]["sha256"],
        "state_sha256"=>checksums["files"]["state.jls"]["sha256"])
end

function csl_solver_settings(cfg, nsweeps, measure_variance)
    s = cfg["solver"]
    return (; seed=s["seed"], initial_linkdim=s["initial_linkdim"], nsweeps,
        maxdim=s["maxdim"], cutoff=0.0, noise=0.0,
        eigsolve_tol=s["eigsolve_tol"], eigsolve_krylovdim=s["eigsolve_krylovdim"],
        eigsolve_maxiter=s["eigsolve_maxiter"], measure_variance)
end

function csl_progress(output, label)
    # Persist only scalar observer fields, never the environment or exception text.
    return function(event)
        if event.kind in (:phase,:sweep)
            progress = csl_dict(event)
            for (key,value) in progress
                value isa Symbol && (progress[key] = string(value))
            end
            progress["worker_elapsed_seconds"] = csl_elapsed()
            progress["unit"] = label
            KagomeDMRG._write_flux_record(joinpath(output,"progress.toml"),progress)
            if event.kind == :sweep
                println("CSL ",label," sweep=",event.sweep," energy=",event.energy,
                    " maxdim=",event.maxlinkdim," elapsed=",round(csl_elapsed();digits=2))
                flush(stdout)
            end
        end
        return nothing
    end
end

function csl_measure!(row, result, phase!, persist!; checkpoint_overlap_error=nothing)
    lattice = result.lattice
    row["energy"] = result.energy
    row["theta"] = result.theta
    row["settings"] = csl_dict(result.settings)
    row["maxlinkdim"] = maxlinkdim(result.psi)
    row["sz_profile"] = copy(result.sz)
    row["column_sz"] = [sum(result.sz[(12x+1):(12x+12)]) for x in 0:2]
    row["sweep_energies"] = copy(result.sweep_energies)
    row["measured_truncation_errors"] = copy(result.max_truncation_errors)
    row["sweep_energy_change"] = max(abs(result.sweep_energies[end]-result.sweep_energies[end-1]),
        abs(result.energy-result.sweep_energies[end]))
    row["state_norm_error"] = abs(inner(result.psi,result.psi)-1)
    row["total_sz_error"] = abs(sum(result.sz))
    row["checkpoint_overlap_measured"] = checkpoint_overlap_error !== nothing
    checkpoint_overlap_error === nothing || (row["checkpoint_overlap_error"] = checkpoint_overlap_error)
    timings = get!(row,"measurement_seconds",Dict{String,Any}())
    phase!(row,"schmidt")
    timings["schmidt"] = @elapsed begin
        schmidt = [schmidt_diagnostics(result.psi,b) for b in CSL_SCHMIDT_BONDS]
        row["schmidt"] = [Dict("bond"=>s.bond,"probabilities"=>s.probabilities,
            "left_q"=>s.left_q,"entropy"=>s.entropy,"mean_left_sz"=>s.mean_left_sz,
            "variance_left_sz"=>s.variance_left_sz) for s in schmidt]
        row["schmidt_density_error"] = maximum(abs(s.mean_left_sz-sum(result.sz[1:s.bond])) for s in schmidt)
        row["schmidt_probability_sum_error"] = maximum(abs(sum(s.probabilities)-1) for s in schmidt)
    end
    persist!()
    phase!(row,"bond_energies")
    timings["bond_energies"] = @elapsed row["bond_energy_profile"] =
        bond_energies(result.psi,lattice,result.theta;gauge=:seam)
    families = bond_families(lattice)
    row["family_energy"] = Dict(string(f)=>sum(row["bond_energy_profile"][families .== f]) for f in (:J1,:J2,:J3))
    row["bond_energy_sum_error"] = abs(sum(row["bond_energy_profile"])-result.energy)
    row["scaled_bond_energy_sum_error"] = row["bond_energy_sum_error"]/max(1,abs(result.energy))
    persist!()
    phase!(row,"scalar_chirality")
    timings["scalar_chirality"] = @elapsed row["chirality_profile"] =
        triangle_chiralities(result.psi,lattice,result.theta;gauge=:seam)
    row["chirality_mean"] = sum(row["chirality_profile"])/length(row["chirality_profile"])
    persist!()
    phase!(row,"variance")
    if result.variance === nothing
        # The saved trial continues to say measure_variance=false. This timed
        # postprocessing diagnostic never rewrites its solver settings/state.
        timings["variance"] = @elapsed begin
            moment = inner(result.H,result.psi,result.H,result.psi)
            row["raw_variance"] = real(moment)-result.energy^2
            row["second_moment_imaginary_part"] = abs(imag(moment))
        end
        row["variance_source"] = "separate_HdaggerH_contraction_after_trial_checkpoint"
        row["second_moment_imaginary_part_measured"] = true
    else
        row["raw_variance"] = result.variance
        row["variance_source"] = "run_dmrg_measured_variance"
        row["second_moment_imaginary_part_measured"] = false
    end
    row["variance_roundoff_scale"] = 100eps(Float64)*max(1.0,result.energy^2)
    failures = [key for (key,limit) in CSL_INTEGRITY_LIMITS
        if !(key=="checkpoint_overlap_error" && !row["checkpoint_overlap_measured"]) &&
            (!isfinite(row[key]) || row[key] > limit)]
    for key in ("sz_profile","bond_energy_profile","chirality_profile","sweep_energies","measured_truncation_errors")
        all(isfinite,row[key]) || push!(failures,"nonfinite_"*key)
    end
    isfinite(row["raw_variance"]) && row["raw_variance"] >= -row["variance_roundoff_scale"] ||
        push!(failures,"invalid_variance")
    if row["second_moment_imaginary_part_measured"]
        row["second_moment_imaginary_part"] <= row["variance_roundoff_scale"] ||
            push!(failures,"second_moment_imaginary_part")
    end
    row["integrity_failures"] = sort(failures)
    row["integrity_passed"] = isempty(failures)
    row["accuracy_status"] = "unestablished"
    row["active_phase"] = "diagnostics_completed"
    persist!()
    return row
end

function csl_save_reload!(output, result, row, phase!, persist!; baseline=result, theta_path=[0.0])
    phase!(row,"trial_checkpoint")
    row["checkpoint_save_seconds"] = @elapsed snapshot = save_checkpoint(output,result;
        baseline,theta_path,status=:trial)
    row["checkpoint"] = relpath(snapshot,output)
    persist!()
    phase!(row,"trial_reload_verification")
    row["checkpoint_reload_seconds"] = @elapsed restored = load_checkpoint(snapshot,result.lattice;
        Q=0,gauge=:seam,status=:trial,expected_theta=result.theta,expected_settings=result.settings)
    row["snapshot"] = csl_snapshot_record(snapshot)
    overlap_error = abs(abs(inner(restored.psi,result.psi))-1)
    row["checkpoint_overlap_error"] = overlap_error
    persist!()
    return snapshot, restored, overlap_error
end

function csl_prepare!(output, cfg, lattice, start, record, phase!, persist!, sources_unchanged!)
    previous = start
    previous_row = nothing
    previous_checkpoint = isempty(cfg["start_checkpoint"]) ? "" : csl_path(cfg["start_checkpoint"])
    start === nothing || start.theta == 0.0 || throw(ArgumentError("prepare requires a zero-flux parent"))
    settings = csl_solver_settings(cfg,cfg["preparation"]["batch_sweeps"],false)
    for batch in 1:cfg["preparation"]["batches"]
        sources_unchanged!()
        row = Dict{String,Any}("batch"=>batch,"status"=>"running","theta"=>0.0,
            "sweeps_in_this_invocation"=>batch*settings.nsweeps,
            "parent_checkpoint"=>isempty(previous_checkpoint) ? "" : csl_record_path(previous_checkpoint))
        push!(record["preparation_batches"],row)
        phase!(row,"DMRG")
        callback = csl_progress(output,"prepare_batch_$(batch)")
        row["dmrg_seconds"] = @elapsed result = previous === nothing ?
            run_dmrg(lattice,0.0;Q=0,settings...,progress_callback=callback) :
            run_dmrg(lattice,0.0;Q=0,sites=previous.sites,psi0=previous.psi,
                settings...,progress_callback=callback)
        snapshot, restored, overlap_error = csl_save_reload!(output,result,row,phase!,persist!)
        # Diagnose the reloaded state against the Hamiltonian used for this batch.
        measured = merge(KagomeDMRG._flux_result(restored),(;H=result.H))
        csl_measure!(row,measured,phase!,persist!;checkpoint_overlap_error=overlap_error)
        if previous !== nothing
            prev = KagomeDMRG._flux_data(previous)
            row["change_from_parent"] = Dict("energy"=>result.energy-prev.energy,
                "overlap"=>abs(inner(previous.psi,result.psi))/(norm(previous.psi)*norm(result.psi)),
                "max_center_sz_change"=>maximum(abs.(result.sz[CSL_CENTER_SITES]-prev.sz[CSL_CENTER_SITES])))
        end
        if previous_row !== nothing
            row["change_from_previous_measured_batch"] = Dict(
                "max_bond_change"=>maximum(abs.(row["bond_energy_profile"]-previous_row["bond_energy_profile"])),
                "max_chirality_change"=>maximum(abs.(row["chirality_profile"]-previous_row["chirality_profile"])),
                "max_entropy_change"=>maximum(abs(row["schmidt"][k]["entropy"]-
                    previous_row["schmidt"][k]["entropy"]) for k in eachindex(CSL_SCHMIDT_BONDS)))
        end
        row["status"] = row["integrity_passed"] ? "completed_accuracy_unestablished" : "integrity_failed"
        record["latest_trial_checkpoint"] = relpath(snapshot,output)
        persist!()
        row["integrity_passed"] || (record["status"]="integrity_failed"; return)
        previous, previous_row, previous_checkpoint = restored, row, snapshot
    end
    record["status"] = "completed_preparation_accuracy_unestablished"
    record["flux_readiness"] = "not_certified_by_preparation_pilot"
    return nothing
end

function csl_continue_solve(saved,lattice,theta;gauge,hz,outputlevel,progress_callback)
    # Exact settings contract of _continuation_solve, with the public progress
    # callback added. continue_flux independently checks these returned settings.
    s = saved.settings
    return run_dmrg(lattice,theta;Q=saved.Q,sites=saved.sites,psi0=saved.psi,
        seed=s.seed,initial_linkdim=s.initial_linkdim,nsweeps=s.nsweeps,
        maxdim=s.maxdim,cutoff=s.cutoff,noise=s.noise,eigsolve_tol=s.eigsolve_tol,
        eigsolve_krylovdim=s.eigsolve_krylovdim,eigsolve_maxiter=s.eigsolve_maxiter,
        measure_variance=s.measure_variance,gauge,hz,outputlevel,progress_callback)
end

function csl_flux!(output,cfg,lattice,start,record,phase!,persist!,sources_unchanged!)
    settings = csl_solver_settings(cfg,cfg["flux"]["sweeps"],true)
    initial_row = Dict{String,Any}("status"=>"running")
    record["flux_start"] = initial_row
    if cfg["start_status"] == "trial"
        start.theta == 0.0 || throw(ArgumentError("trial flux parent must be at zero flux"))
        # A new optimization and new variance-enabled snapshot establish a new
        # measured baseline. Never label an old prepare snapshot as eligible.
        phase!(initial_row,"new_variance_enabled_zero_flux_batch")
        initial_row["dmrg_seconds"] = @elapsed initial = run_dmrg(lattice,0.0;
            Q=0,sites=start.sites,psi0=start.psi,settings...,
            progress_callback=csl_progress(output,"flux_initial_zero"))
        _, restored, reload_error = csl_save_reload!(output,initial,initial_row,phase!,persist!)
        csl_measure!(initial_row,KagomeDMRG._flux_result(restored),phase!,persist!;
            checkpoint_overlap_error=reload_error)
        initial_path = [0.0]
        continuation_start = initial
        initial_row["baseline_origin"] = "new_zero_flux_variance_enabled_batch"
    else
        # Restarting an accepted trajectory retains its entire measured baseline
        # and settings. A config mismatch requires a separate new experiment.
        actual = start.settings
        for (key,value) in pairs(settings)
            expected = key == :maxdim ? [value] : value
            getproperty(actual,key) == expected || throw(ArgumentError("accepted restart solver mismatch"))
        end
        actual.measure_variance && actual.nsweeps >= 2 || throw(ArgumentError("ineligible accepted start"))
        csl_measure!(initial_row,KagomeDMRG._flux_result(start),phase!,persist!)
        initial_path = copy(start.theta_path)
        continuation_start = csl_path(cfg["start_checkpoint"])
        initial_row["baseline_origin"] = "original_measured_baseline_from_accepted_checkpoint"
    end
    initial_row["status"] = initial_row["integrity_passed"] ? "awaiting_core_policy" : "integrity_failed"
    initial_row["theta_path"] = initial_path
    persist!()
    initial_row["integrity_passed"] || (record["status"]="integrity_failed"; return)
    triangles = oriented_triangles(lattice)
    center_bonds = [k for (k,b) in enumerate(lattice.bonds) if
        b.i in CSL_CENTER_SITES && b.j in CSL_CENTER_SITES]
    center_triangles = [k for (k,t) in enumerate(triangles) if all(i->i in CSL_CENTER_SITES,t.sites)]
    (isempty(center_bonds) || isempty(center_triangles)) && error("empty central diagnostics")
    record["extra_branch_monitor"] = Dict("bond_indices"=>center_bonds,
        "triangle_indices"=>center_triangles,"site_selection"=>CSL_CENTER_SITES,
        "meaning"=>"central_column_monitor_not_a_demonstrated_bulk",
        "max_bond_change"=>cfg["flux"]["max_bond_change"],
        "max_chirality_change"=>cfg["flux"]["max_chirality_change"],
        "expected_pump_used"=>false)
    # Only keys supplied by the core's last accepted checkpoint are ever read.
    # Observed/rejected candidates may be cached but can never seed a solve.
    cache = Dict{Tuple,Any}(Tuple(initial_path)=>initial_row)
    function solver(saved,lat,theta;gauge,hz,outputlevel)
        sources_unchanged!()
        parent_key = Tuple(saved.theta_path)
        haskey(cache,parent_key) || throw(ArgumentError("missing accepted-path diagnostics"))
        row = Dict{String,Any}("status"=>"running","theta"=>theta,
            "from_theta"=>saved.theta,"parent_theta_path"=>copy(saved.theta_path),
            "theta_path"=>[saved.theta_path;theta],"core_acceptance"=>"not_assigned_by_worker")
        push!(record["flux_observations"],row)
        try
            phase!(row,"DMRG")
            row["dmrg_seconds"] = @elapsed result = csl_continue_solve(saved,lat,theta;
                gauge,hz,outputlevel,progress_callback=csl_progress(output,"flux_trial_$(length(record["flux_observations"]))"))
            _, restored, reload_error = csl_save_reload!(joinpath(output,"extra_observer_trials"),
                result,row,phase!,persist!;baseline=saved.baseline,theta_path=[saved.theta_path;theta])
            # The snapshot path in this nested staging area is ROOT-relative in
            # row.snapshot.path. It remains a trial even if core later accepts it.
            row["checkpoint_root"] = "extra_observer_trials"
            csl_measure!(row,KagomeDMRG._flux_result(restored),phase!,persist!;
                checkpoint_overlap_error=reload_error)
            row["integrity_passed"] || throw(ArgumentError("extra diagnostic integrity failure"))
            previous = cache[parent_key]
            row["max_center_bond_change"] = maximum(abs.(row["bond_energy_profile"][center_bonds]-
                previous["bond_energy_profile"][center_bonds]))
            row["max_center_chirality_change"] = maximum(abs.(row["chirality_profile"][center_triangles]-
                previous["chirality_profile"][center_triangles]))
            reasons = String[]
            row["max_center_bond_change"] <= cfg["flux"]["max_bond_change"] || push!(reasons,"center_bond_change")
            row["max_center_chirality_change"] <= cfg["flux"]["max_chirality_change"] || push!(reasons,"center_chirality_change")
            row["extra_rejection_reasons"] = reasons
            row["status"] = isempty(reasons) ? "returned_to_core_policy" : "extra_branch_rejected"
            cache[Tuple([saved.theta_path;theta])] = row
            persist!()
            # ErrorException is the existing core's retry category: it reloads last
            # accepted and halves the step. Detailed extra reasons remain above.
            isempty(reasons) || throw(ErrorException("extra_branch_diagnostic_rejection"))
            return result
        catch exception
            row["status"] == "running" && (row["status"] = "solver_or_diagnostic_failed")
            row["exception_type"] = string(nameof(typeof(exception)))
            persist!()
            rethrow()
        end
    end
    persist!()
    p = cfg["policy"]
    policy = FluxPolicy(; (Symbol(k)=>v for (k,v) in p)...)
    f = cfg["flux"]
    result = continue_flux(lattice,f["targets"];start=continuation_start,
        output_root=joinpath(output,"continuation"),policy,initial_step=f["initial_step"],
        min_step=f["min_step"],max_trials=f["max_trials"],cuts=CSL_CUTS,
        diagnostic_bonds=CSL_SCHMIDT_BONDS,bulk_sites=CSL_CENTER_SITES,
        gauge=:seam,Q=0,outputlevel=0,_solver=solver)
    record["continuation"] = Dict("status"=>string(result.status),"reason"=>result.reason,
        "theta_path"=>result.theta_path,"journal"=>relpath(result.output_path,output),
        "last_accepted_checkpoint"=>result.last_checkpoint === nothing ? "" : relpath(result.last_checkpoint,output),
        "accepted_checkpoints"=>relpath.(result.accepted_checkpoints,Ref(output)))
    record["status"] = result.status == :completed ? "completed_finite_flux_limits_only" : "unresolved_flux"
    record["flux_start"]["status"] = "see_core_initial_policy_decision"
    return nothing
end

function csl_main(output,config_path)
    output, config_path = abspath(output), abspath(config_path)
    isdir(output) && readdir(output) == ["execution.toml"] ||
        throw(ArgumentError("use run_research.py and a fresh output directory"))
    execution = TOML.parsefile(joinpath(output,"execution.toml"))
    get(execution,"status","") == "running" || throw(ArgumentError("launcher is not running"))
    wall = get(execution,"wall_limit_seconds",Inf)
    wall isa Real && isfinite(wall) && 0 < wall <= 600 || throw(ArgumentError("wall limit must be at most 600 s"))
    Threads.nthreads() == BLAS.get_num_threads() == 1 || throw(ArgumentError("one Julia/BLAS thread required"))
    cfg_sha = csl_sha(config_path)
    cfg = csl_read_config(config_path)
    source_hashes = csl_extra_hashes()
    sources_unchanged!() = begin
        csl_sha(config_path) == cfg_sha && csl_extra_hashes() == source_hashes ||
            throw(ArgumentError("worker configuration or source changed"))
        KagomeDMRG._execution_identity()
        nothing
    end
    lattice = kagome_j1j2j3_cylinder(3,4;CSL_COUPLINGS...)
    record = Dict{String,Any}("schema_version"=>1,"status"=>"running",
        "recorded_at_utc"=>string(now(UTC)),"scope"=>"known_CSL_model_finite_measurement_control",
        "mode"=>cfg["mode"],"config"=>cfg,"config_sha256"=>cfg_sha,
        "N"=>36,"Q"=>0,"couplings"=>csl_dict(CSL_COUPLINGS),
        "configuration"=>KagomeDMRG._checkpoint_configuration(lattice,:seam,nothing;Q=0),
        "code"=>KagomeDMRG._checkpoint_provenance(),"extra_source_sha256"=>source_hashes,
        "runtime"=>merge(KagomeDMRG._checkpoint_runtime(),Dict("julia_threads"=>1,"blas_threads"=>1)),
        "wall_limit_seconds"=>wall,"execution_record"=>"execution.toml",
        "time_limit_enforcement"=>"external_launcher_including_startup_and_single_solves",
        "cuts"=>CSL_CUTS,"schmidt_bonds"=>CSL_SCHMIDT_BONDS,"center_sites"=>CSL_CENTER_SITES,
        "bulk_status"=>"three_columns_do_not_establish_a_bulk",
        "accuracy_status"=>"unestablished","phase_identification"=>"not_attempted",
        "quantized_pump"=>"not_established","expected_pump_used_for_acceptance"=>false,
        "initial_state_convention"=>"complex_fixed_Q_random_seed_or_explicit_saved_parent",
        "added_chirality_term"=>false,"uniform_zeeman_field"=>0.0,
        "integrity_limits"=>CSL_INTEGRITY_LIMITS,
        "triangles"=>[Dict("sites"=>collect(t.sites),"image_y"=>collect(t.image_y),
            "kind"=>string(t.kind)) for t in oriented_triangles(lattice)],
        "literature"=>Dict("paper"=>"https://www.nature.com/articles/srep06317",
            "author_manuscript"=>"https://arxiv.org/abs/1312.4519",
            "comparison"=>"same_couplings_and_hexagon_opposite_J3_not_same_size_or_convergence",
            "published_flux_case"=>"N288_Lx24_Ly4_Jprime0.5_complex_U1_up_to_5000_states"),
        "preparation_batches"=>Dict{String,Any}[],"flux_observations"=>Dict{String,Any}[])
    function persist!()
        record["worker_elapsed_seconds"] = csl_elapsed()
        KagomeDMRG._write_flux_record(joinpath(output,"validation.toml"),record)
    end
    function phase!(row,phase)
        row["active_phase"] = phase
        persist!()
        println("CSL phase=",phase," elapsed=",round(csl_elapsed();digits=2))
        flush(stdout)
    end
    persist!()
    for relative in CSL_EXTRA_SOURCES
        destination = joinpath(output,"analysis-sources",relative)
        mkpath(dirname(destination))
        cp(joinpath(CSL_ROOT,relative),destination)
    end
    cp(config_path,joinpath(output,"analysis-sources","config.toml"))
    parent_files = Dict{String,String}()
    try
        start = nothing
        if !isempty(cfg["start_checkpoint"])
            path = csl_path(cfg["start_checkpoint"])
            start = load_checkpoint(path,lattice;Q=0,gauge=:seam,status=Symbol(cfg["start_status"]))
            parent_files = Dict(joinpath(path,name)=>csl_sha(joinpath(path,name))
                for name in ("metadata.toml","state.jls","checksums.toml"))
            record["parent"] = merge(csl_snapshot_record(path),Dict("status"=>cfg["start_status"],
                "theta_path"=>copy(start.theta_path),"solver"=>csl_dict(start.settings),
                "energy"=>start.diagnostics.energy,"variance_measured"=>start.settings.measure_variance))
            persist!()
        end
        if cfg["mode"] == "prepare"
            csl_prepare!(output,cfg,lattice,start,record,phase!,persist!,sources_unchanged!)
        else
            csl_flux!(output,cfg,lattice,start,record,phase!,persist!,sources_unchanged!)
        end
    catch exception
        record["status"] = exception isa InterruptException ? "interrupted" : "failed"
        record["exception_type"] = string(nameof(typeof(exception)))
        # Arbitrary exception messages can contain host/path information.
        persist!()
        exception isa InterruptException && rethrow()
        exception isa OutOfMemoryError && rethrow()
    finally
        try
            sources_unchanged!()
            record["configuration_and_source_unchanged"] = true
        catch exception
            record["configuration_and_source_unchanged"] = false
            record["status"] = "failed_source_or_configuration_changed"
            record["source_check_exception_type"] = string(nameof(typeof(exception)))
        end
        record["parent_unchanged"] = all(isfile(path) && csl_sha(path)==hash for (path,hash) in parent_files)
        record["parent_unchanged"] || (record["status"]="failed_parent_changed")
        persist!()
    end
    println("CSL worker status=",record["status"])
    flush(stdout)
    return record
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 2 || error("usage: research_csl.jl OUTPUT CONFIG.toml")
    csl_result = csl_main(ARGS[1],ARGS[2])
    exit(startswith(csl_result["status"],"completed_") ? 0 : 2)
end
