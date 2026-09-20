#!/usr/bin/env julia
# One declared NN case. Keep backend identity distinct from this analysis driver.
const WORKER_START = time_ns()
using KagomeDMRG, ITensors, ITensorMPS, LinearAlgebra, SHA, TOML, Dates
BLAS.set_num_threads(1)
ITensors.disable_threaded_blocksparse()
const ROOT = normpath(joinpath(@__DIR__,".."))
const BASE_SOURCES = ("examples/research_static.jl", "examples/run_research.py",
    "examples/run_static27.py", "examples/run_static18.py")
const INTEGRITY_LIMITS = Dict("norm"=>1e-12,"charge"=>1e-12,"zz_diagonal"=>1e-10,
    "pm_diagonal"=>1e-10,"zz_symmetry"=>1e-10,"zz_imaginary"=>1e-10,
    "pm_hermiticity"=>1e-10,"fixed_charge_zz"=>1e-10,"bond_sum"=>1e-9,
    "schmidt_density"=>1e-10,"schmidt_variance"=>1e-9)
elapsed() = (time_ns()-WORKER_START)/1e9
file_hash(path) = bytes2hex(sha256(read(path)))
rootpath(path) = normpath(isabspath(path) ? path : joinpath(ROOT,path))
asdict(value) = Dict(string(k)=>v for (k,v) in pairs(value))
real_rows(M) = [real.(collect(row)) for row in eachrow(M)]
imag_rows(M) = [imag.(collect(row)) for row in eachrow(M)]

function read_config(path)
    cfg = TOML.parsefile(path)
    allowed = Set(["schema_version","case_id","Lx","Ly","Q","seed","maxdim",
        "batches","batch_sweeps","initialization","parent_checkpoint","parent_record",
        "parent_record_sha256","parent_metadata_sha256","initial_completed_sweeps",
        "backend_snapshot_origin","stationarity","helper"])
    all(k in allowed for k in keys(cfg)) || error("unknown config key")
    cfg["schema_version"]==1 && 1<=cfg["batches"]<=4 && cfg["batch_sweeps"] in (1,2) &&
        3<=cfg["Lx"]<=6 && cfg["Ly"]==3 && cfg["maxdim"] in (128,256,512,1024) || error("unsupported bounded case")
    target_sector(3cfg["Lx"]*cfg["Ly"];Q=cfg["Q"])
    cfg["initialization"] in ("random","resume","period9","period27") || error("unknown initialization")
    initial_sweeps = get(cfg,"initial_completed_sweeps",0)
    initial_sweeps isa Integer && !(initial_sweeps isa Bool) && initial_sweeps>=0 ||
        error("initial_completed_sweeps must be a nonnegative integer")
    cfg["initialization"]=="resume" || initial_sweeps==0 || error("fresh run cannot have prior sweeps")
    limits = cfg["stationarity"]
    Set(keys(limits))==Set(["energy_per_site","sz_profile","bond_profile","variance_per_site","truncation"]) &&
        all(v isa Real && isfinite(v) && v>0 for v in values(limits)) || error("invalid numerical comparison limits")
    return cfg
end

function parent_state(cfg,lattice,record)
    path = rootpath(cfg["parent_checkpoint"])
    evidence_path = rootpath(cfg["parent_record"])
    file_hash(evidence_path)==cfg["parent_record_sha256"] || error("parent evidence hash mismatch")
    file_hash(joinpath(path,"metadata.toml"))==cfg["parent_metadata_sha256"] || error("parent metadata hash mismatch")
    evidence = TOML.parsefile(evidence_path)
    previous = if haskey(evidence,"batches")
        only(filter(r->haskey(r,"checkpoint") && rootpath(r["checkpoint"])==path,evidence["batches"]))
    else
        evidence
    end
    old_sweeps = get(previous,"cumulative_sweeps",get(previous,"completed_sweeps",-1))
    old_sweeps==cfg["initial_completed_sweeps"] && rootpath(previous["checkpoint"])==path || error("parent sweep/path mismatch")
    get(previous,"integrity_passed",get(previous,"state_integrity_passed",false)) || error("parent integrity not established")
    snapshot_files = Dict(relpath(joinpath(path,f),ROOT)=>file_hash(joinpath(path,f))
        for f in ("metadata.toml","state.jls","checksums.toml"))
    snapshot_files[relpath(evidence_path,ROOT)] = cfg["parent_record_sha256"]
    if haskey(evidence,"validation_record")
        validation_path = normpath(joinpath(dirname(evidence_path),evidence["validation_record"]))
        validation_hash = evidence["validation_record_sha256"]
        file_hash(validation_path)==validation_hash || error("parent validation hash mismatch")
        snapshot_files[relpath(validation_path,ROOT)] = validation_hash
        validation = TOML.parsefile(validation_path)
        validation["N"]==nsites(lattice) && validation["Q"]==cfg["Q"] || error("wrong parent validation sector")
        previous = merge(previous,Dict("sz_profile"=>validation["state"]["sz_profile"],
            "bond_energy"=>validation["bond_observables"]["energy_profile"]))
    end
    meta = TOML.parsefile(joinpath(path,"metadata.toml"))
    expected = KagomeDMRG._checkpoint_settings(meta["settings"])
    saved = load_checkpoint(path,lattice;Q=cfg["Q"],status=:trial,expected_theta=0.0,expected_settings=expected)
    record["parent"] = Dict("checkpoint"=>relpath(path,ROOT),"record"=>relpath(evidence_path,ROOT),
        "record_sha256"=>cfg["parent_record_sha256"],"checkpoint_sha256"=>snapshot_files,
        "initial_completed_sweeps"=>old_sweeps,"settings"=>asdict(saved.settings),
        "loaded_energy"=>saved.diagnostics.energy,"loaded_norm_error"=>abs(norm(saved.psi)-1),
        "strict_source_runtime_configuration_load"=>true)
    return saved,previous,snapshot_files
end

function measure_state(result,lattice,record,row,persist)
    N,Q = nsites(lattice),result.Q
    start = elapsed()
    row["active_phase"] = "correlations"
    persist()
    corr = spin_correlations(result.psi)
    bonds = [b.Jz*real(corr.zz[b.i,b.j])+b.Jxy*real(corr.pm[b.i,b.j]) for b in lattice.bonds]
    row["bond_energy"] = bonds
    row["correlations"] = Dict("zz_real"=>real_rows(corr.zz),"zz_imag"=>imag_rows(corr.zz),
        "pm_real"=>real_rows(corr.pm),"pm_imag"=>imag_rows(corr.pm))
    row["correlation_seconds"] = elapsed()-start
    row["active_phase"] = "variance"
    persist()
    variance_seconds = @elapsed second_moment = inner(result.H,result.psi,result.H,result.psi)
    variance = real(second_moment)-result.energy^2
    row["HdaggerH_real"] = real(second_moment)
    row["HdaggerH_imag"] = imag(second_moment)
    row["variance"] = variance
    row["variance_per_site"] = variance/N
    row["variance_seconds"] = variance_seconds
    row["variance_roundoff_scale"] = 100eps(Float64)*max(1,result.energy^2)
    row["active_phase"] = "Schmidt"
    persist()
    schmidt = Dict{String,Any}[]
    for cut in 1:lattice.Lx-1
        b = 3lattice.Ly*cut
        d = schmidt_diagnostics(result.psi,b)
        push!(schmidt,merge(asdict(d),Dict("cut"=>cut,"bond"=>b,
            "density_error"=>abs(d.mean_left_sz-sum(result.sz[1:b])),
            "variance_error"=>abs(d.variance_left_sz-(real(sum(corr.zz[1:b,1:b]))-sum(result.sz[1:b])^2)))))
    end
    row["schmidt"] = schmidt
    errors = Dict("norm"=>abs(norm(result.psi)-1),"charge"=>abs(sum(result.sz)-Q/2),
        "zz_diagonal"=>maximum(abs.(diag(corr.zz).-0.25)),
        "pm_diagonal"=>maximum(abs.(diag(corr.pm).-(0.5 .+result.sz))),
        "zz_symmetry"=>maximum(abs.(corr.zz-transpose(corr.zz))),
        "zz_imaginary"=>maximum(abs.(imag.(corr.zz))),
        "pm_hermiticity"=>maximum(abs.(corr.pm-corr.pm')),
        "fixed_charge_zz"=>maximum(abs.(vec(sum(corr.zz;dims=2))-(Q/2).*result.sz)),
        "bond_sum"=>abs(sum(bonds)-result.energy),
        "schmidt_density"=>maximum(r["density_error"] for r in schmidt),
        "schmidt_variance"=>maximum(r["variance_error"] for r in schmidt))
    row["integrity_errors"] = errors
    row["integrity_passed"] = all(isfinite(v) && v<=INTEGRITY_LIMITS[k] for (k,v) in errors) &&
        isfinite(variance) && variance>=-row["variance_roundoff_scale"] &&
        isfinite(second_moment) && abs(imag(second_moment))<=row["variance_roundoff_scale"] &&
        all(isfinite,result.sz) && all(isfinite,bonds) && all(isfinite,corr.zz) && all(isfinite,corr.pm) &&
        all(x->isfinite(x) && 0<=x<=1,result.max_truncation_errors) && row["checkpoint_verified"]
    row["diagnostic_seconds"] = elapsed()-start
    row["status"] = row["integrity_passed"] ? "diagnostics_passed_accuracy_separate" : "failed_integrity"
    row["active_phase"] = "completed"
    return nothing
end

function main(output,config_path)
    isdir(output) && readdir(output)==["execution.toml"] || error("use fresh launcher output")
    execution = TOML.parsefile(joinpath(output,"execution.toml"))
    execution["status"]=="running" && 0<execution["wall_limit_seconds"]<=600 || error("invalid supervisor")
    Threads.nthreads()==BLAS.get_num_threads()==1 || error("one Julia/BLAS thread required")
    cfg = read_config(config_path)
    sources = collect(BASE_SOURCES)
    if cfg["initialization"] in ("period9","period27")
        push!(sources,"examples/order_seed_tools.jl")
    end
    before = Dict(p=>file_hash(joinpath(ROOT,p)) for p in sources)
    config_hash = file_hash(config_path)
    lattice = kagome_cylinder(cfg["Lx"],cfg["Ly"])
    N,Q = nsites(lattice),cfg["Q"]
    settings = (; seed=cfg["seed"],nsweeps=cfg["batch_sweeps"],maxdim=cfg["maxdim"],
        cutoff=0.0,noise=0.0,eigsolve_tol=1e-11,eigsolve_krylovdim=40,eigsolve_maxiter=20,measure_variance=false)
    record = Dict{String,Any}("schema_version"=>1,"status"=>"running","case_id"=>cfg["case_id"],
        "recorded_at_utc"=>string(now(UTC)),"N"=>N,"Q"=>Q,"theta"=>0.0,"configuration"=>KagomeDMRG._checkpoint_configuration(lattice,:seam,nothing;Q),
        "config"=>cfg,"config_path"=>relpath(config_path,ROOT),"config_sha256"=>config_hash,
        "analysis_source_sha256"=>before,"code"=>KagomeDMRG._checkpoint_provenance(),
        "backend_project"=>relpath(pkgdir(KagomeDMRG),ROOT),
        "backend_snapshot_origin"=>get(cfg,"backend_snapshot_origin","current_workspace"),
        "git_revision_semantics"=>"execution_directory_git_lookup_not_a_replacement_for_snapshot_hashes",
        "runtime"=>merge(KagomeDMRG._checkpoint_runtime(),Dict("julia_threads"=>1,"blas_threads"=>1)),
        "solver"=>asdict(settings),"integrity_limits"=>INTEGRITY_LIMITS,
        "wall_limit_seconds"=>execution["wall_limit_seconds"],
        "new_ED"=>false,"phase_identification"=>"not_attempted","bulk_plateau"=>"not_established",
        "variance_interpretation"=>"H_dispersion_not_ground_energy_error_bound",
        "stationarity_interpretation"=>"fixed_chi_observation_requires_chi_seed_size_checks",
        "chirality"=>"not_measured_in_this_NN_worker","batches"=>Dict{String,Any}[])
    persist() = KagomeDMRG._write_flux_record(joinpath(output,"validation.toml"),record)
    persist()
    # Preserve exact research code even if a later campaign run revises it.
    for p in sources
        destination = joinpath(output,"analysis-sources",p)
        mkpath(dirname(destination))
        cp(joinpath(ROOT,p),destination)
    end
    cp(config_path,joinpath(output,"analysis-sources","config.toml"))
    parent_hashes = Dict{String,String}()
    completed_sweeps = get(cfg,"initial_completed_sweeps",0)
    sites,psi0,previous = spin_sites(lattice),nothing,nothing
    previous_maxdim = nothing
    try
        if cfg["initialization"]=="resume"
            saved,previous,parent_hashes = parent_state(cfg,lattice,record)
            sites,psi0 = saved.sites,saved.psi
            previous_maxdim = length(saved.settings.maxdim)==1 ? only(saved.settings.maxdim) : nothing
        elseif cfg["initialization"] in ("period9","period27")
            seed_state = Base.invokelatest(make_order_seed,lattice,sites;kind=Symbol(cfg["initialization"]),Q,seed=cfg["seed"])
            psi0 = seed_state.psi
            record["seed_metadata"] = seed_state.metadata
        end
        persist()
        for batch in 1:cfg["batches"]
            row = Dict{String,Any}("batch"=>batch,"status"=>"running","active_phase"=>"DMRG",
                "cumulative_sweeps"=>completed_sweeps+settings.nsweeps)
            push!(record["batches"],row)
            function progress(event)
                if event.kind==:sweep || event.kind==:phase
                    row["last_event"] = Dict(string(k)=>(v isa Symbol ? string(v) : v) for (k,v) in pairs(event))
                    record["worker_elapsed_seconds"] = elapsed()
                    persist()
                    if event.kind==:sweep
                        println(cfg["case_id"]," sweep=",completed_sweeps+event.sweep," E=",event.energy," chi=",event.maxlinkdim)
                        flush(stdout)
                    end
                end
            end
            dmrg_seconds = @elapsed result = run_dmrg(lattice,0.0;Q,sites,psi0,settings...,progress_callback=progress)
            row["dmrg_seconds"] = dmrg_seconds
            row["energy"] = result.energy
            row["energy_per_site"] = result.energy/N
            row["sz_profile"] = result.sz
            row["column_sz"] = [sum(result.sz[3lattice.Ly*x+1:3lattice.Ly*(x+1)]) for x in 0:lattice.Lx-1]
            row["sweep_energies"] = result.sweep_energies
            row["measured_truncation_errors"] = result.max_truncation_errors
            row["maxlinkdim"] = maxlinkdim(result.psi)
            row["last_optimizer_energy_error"] = abs(result.energy-last(result.sweep_energies))
            row["within_batch_sweep_change_measured"] = length(result.sweep_energies)>=2
            if row["within_batch_sweep_change_measured"]
                row["last_sweep_energy_change"] = max(abs(diff(result.sweep_energies)[end]),row["last_optimizer_energy_error"])
            end
            row["active_phase"] = "checkpoint"
            persist()
            snapshot = save_checkpoint(output,result;baseline=result,theta_path=[0.0],status=:trial)
            row["checkpoint"] = relpath(snapshot,ROOT)
            row["checkpoint_metadata_sha256"] = file_hash(joinpath(snapshot,"metadata.toml"))
            row["checkpoint_payload_sha256"] = file_hash(joinpath(snapshot,"state.jls"))
            persist()
            saved = load_checkpoint(snapshot,lattice;Q,status=:trial,expected_theta=0.0,expected_settings=result.settings)
            row["checkpoint_verified"] = abs(abs(inner(saved.psi,result.psi))-1)<=1e-12
            row["status"] = "checkpoint_saved_diagnostics_pending"
            persist()
            measure_state(result,lattice,record,row,persist)
            if previous!==nothing
                row["same_chi_comparison"] = previous_maxdim==settings.maxdim
                row["comparison_kind"] = row["same_chi_comparison"] ? "same_chi" : "chi_change"
                row["energy_change_from_previous"] = result.energy-previous["energy"]
                if haskey(previous,"sz_profile") && haskey(previous,"bond_energy")
                    row["max_sz_change"] = maximum(abs.(result.sz-previous["sz_profile"]))
                    row["max_bond_change"] = maximum(abs.(row["bond_energy"]-previous["bond_energy"]))
                    limits = cfg["stationarity"]
                    row["stationarity_evaluated"] = row["same_chi_comparison"]
                    row["stationarity_passed"] = row["same_chi_comparison"] && max(abs(row["energy_change_from_previous"]),
                        get(row,"last_sweep_energy_change",row["last_optimizer_energy_error"]))/N<=limits["energy_per_site"] &&
                        row["max_sz_change"]<=limits["sz_profile"] && row["max_bond_change"]<=limits["bond_profile"] &&
                        abs(row["variance_per_site"])<=limits["variance_per_site"] &&
                        last(row["measured_truncation_errors"])<=limits["truncation"]
                end
            end
            row["integrity_passed"] || error("observable consistency failed")
            completed_sweeps += settings.nsweeps
            previous,sites,psi0 = row,result.sites,result.psi
            previous_maxdim = settings.maxdim
            persist()
        end
        record["status"] = "completed_bounded_case_accuracy_separate"
    catch exception
        record["status"] = "failed_or_interrupted"
        record["exception_type"] = string(typeof(exception))
        rethrow()
    finally
        record["worker_elapsed_seconds"] = elapsed()
        record["completed_sweeps"] = completed_sweeps
        record["sources_unchanged"] = all(file_hash(joinpath(ROOT,p))==h for (p,h) in before) && file_hash(config_path)==config_hash
        try
            record["backend_unchanged"] = KagomeDMRG._checkpoint_provenance()["source_sha256"]==record["code"]["source_sha256"]
        catch exception
            record["backend_unchanged"] = false
            record["backend_check_exception_type"] = string(typeof(exception))
        end
        record["parent_unchanged"] = all(file_hash(joinpath(ROOT,p))==h for (p,h) in parent_hashes)
        record["sources_unchanged"] && record["backend_unchanged"] && record["parent_unchanged"] || (record["status"]="source_or_parent_changed")
        persist()
    end
    return record
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==2 || error("expected OUTPUT CONFIG.toml")
    cfg = read_config(abspath(ARGS[2]))
    if cfg["initialization"] in ("period9","period27")
        include("order_seed_tools.jl")
    end
    record = main(abspath(ARGS[1]),abspath(ARGS[2]))
    println(record["status"])
    record["status"]=="completed_bounded_case_accuracy_separate" || exit(1)
end
