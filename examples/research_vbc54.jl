#!/usr/bin/env julia
# Reuse the existing observable/integrity definitions without launching its driver.
include("research_static.jl")
include("vbc_motif_tools.jl")

const VBC_SOURCES = ("examples/research_vbc54.jl", "examples/vbc_motif_tools.jl",
    "examples/research_static.jl", "examples/run_vbc54.py", "examples/run_static27.py",
    "examples/run_static18.py")
const VBC_PRECISION = Dict("energy_per_site"=>1e-6,"sz_profile"=>1e-4,
    "bond_profile"=>1e-4,"variance_per_site"=>1e-5,"truncation"=>1e-6)

function vbc_config(path)
    cfg = TOML.parsefile(path)
    Set(keys(cfg))==Set(["schema_version","case_id","Lx","Ly","Q","seed",
        "origin","branches","wall_seconds","question","stages","precision"]) || error("unknown configuration keys")
    cfg["schema_version"]==1 && (cfg["Lx"],cfg["Ly"],cfg["Q"])==(6,3,6) || error("expected N54 Q6")
    cfg["branches"]==["random","hourglass","windmill"] && cfg["origin"]==[0,0] || error("unsupported preparation comparison")
    cfg["wall_seconds"]==2700 && cfg["precision"]==VBC_PRECISION || error("changed budget or precision criteria")
    expected = [("prepare",0.3,[32,64],false),("reduce",0.1,[128,128],false),
        ("release",0.0,[128,128],true),("relax",0.0,[128,128],true)]
    length(cfg["stages"])==4 || error("expected four declared stages")
    for (s,(name,lambda,dims,measure)) in zip(cfg["stages"],expected)
        Set(keys(s))==Set(["name","lambda","maxdim","measure"]) &&
            (s["name"],s["lambda"],s["maxdim"],s["measure"])==(name,lambda,dims,measure) || error("changed stage protocol")
    end
    return cfg
end

function vbc_final_model(lattice,canonical)
    lattice.sites==canonical.sites && length(lattice.bonds)==length(canonical.bonds) || return false
    return all((a.i,a.j,a.wy,a.Jxy,a.Jz)==(b.i,b.j,b.wy,1.0,1.0)
        for (a,b) in zip(lattice.bonds,canonical.bonds))
end

function vbc_compare_precision!(row,previous,N)
    row["energy_change_from_previous"] = row["energy"]-previous["energy"]
    row["max_sz_change"] = maximum(abs.(row["sz_profile"]-previous["sz_profile"]))
    row["max_bond_change"] = maximum(abs.(row["bond_energy"]-previous["bond_energy"]))
    values = Dict("energy_per_site"=>max(abs(row["energy_change_from_previous"]),
        row["last_sweep_energy_change"],row["last_optimizer_energy_error"])/N,
        "sz_profile"=>row["max_sz_change"],"bond_profile"=>row["max_bond_change"],
        "variance_per_site"=>abs(row["variance_per_site"]),
        "truncation"=>last(row["measured_truncation_errors"]))
    row["precision_values"] = values
    row["precision_passed"] = Dict(k=>isfinite(v) && v<=VBC_PRECISION[k] for (k,v) in values)
    row["all_precision_conditions_passed"] = all(Base.values(row["precision_passed"]))
    row["comparison_kind"] = "same_chi_same_unpinned_H_additional_two_sweeps"
end

function vbc_main(output,config_path)
    isdir(output) && readdir(output)==["execution.toml"] || error("use a fresh supervised output")
    cfg = vbc_config(config_path)
    execution = TOML.parsefile(joinpath(output,"execution.toml"))
    execution["status"]=="running" && execution["wall_limit_seconds"]==cfg["wall_seconds"] || error("invalid supervisor")
    Threads.nthreads()==BLAS.get_num_threads()==1 || error("one Julia/BLAS thread required")
    lattice = kagome_cylinder(6,3)
    sites = spin_sites(lattice)
    initial = initial_mps(sites;Q=6,seed=cfg["seed"],linkdim=4)
    sources = Dict(p=>file_hash(joinpath(ROOT,p)) for p in VBC_SOURCES)
    config_hash = file_hash(config_path)
    record = Dict{String,Any}("schema_version"=>1,"case_id"=>cfg["case_id"],
        "status"=>"running","recorded_at_utc"=>string(now(UTC)),"N"=>54,"Q"=>6,"theta"=>0.0,
        "config"=>cfg,"config_sha256"=>config_hash,"analysis_source_sha256"=>sources,
        "configuration"=>KagomeDMRG._checkpoint_configuration(lattice,:seam,nothing;Q=6),
        "code"=>KagomeDMRG._checkpoint_provenance(),
        "runtime"=>merge(KagomeDMRG._checkpoint_runtime(),Dict("julia_threads"=>1,"blas_threads"=>1)),
        "wall_limit_seconds"=>cfg["wall_seconds"],"integrity_limits"=>INTEGRITY_LIMITS,
        "initial_state"=>Dict("kind"=>"same_complex_random_fixed_charge_mps_for_all_branches",
            "seed"=>cfg["seed"],"linkdim"=>4,"norm_error"=>abs(norm(initial)-1),
            "sz_profile"=>sz_profile(initial)),
        "comparison_scope"=>"finite_chi_preparation_sensitivity_not_phase_identification",
        "preparation_rule"=>"Jxy=Jz=1+lambda*weight; random branch has lambda=0 at every stage",
        "stage_order"=>"round_robin_across_branches_before_advancing_stage",
        "variance_interpretation"=>"H_dispersion_not_ground_energy_error_bound",
        "templates"=>Dict{String,Any}(),"batches"=>Dict{String,Any}[])
    persist() = KagomeDMRG._write_flux_record(joinpath(output,"validation.toml"),record)
    persist()
    for p in VBC_SOURCES
        dest = joinpath(output,"analysis-sources",p)
        mkpath(dirname(dest)); cp(joinpath(ROOT,p),dest)
    end
    cp(config_path,joinpath(output,"analysis-sources/config.toml"))
    states = Dict(branch=>deepcopy(initial) for branch in cfg["branches"])
    previous = Dict{String,Any}()
    templates = Dict{String,Any}()
    for branch in ("hourglass","windmill")
        template = vbc_bond_template(lattice;kind=Symbol(branch),origin=Tuple(cfg["origin"]))
        templates[branch] = template
        record["templates"][branch] = template.metadata
        record["templates"][branch]["weights"] = template.weights
    end
    persist()
    try
        for (stage_index,stage) in enumerate(cfg["stages"]), branch in cfg["branches"]
            lambda = branch=="random" ? 0.0 : stage["lambda"]
            model = branch=="random" ? lattice : vbc_prepared_lattice(lattice,templates[branch],lambda)
            unpinned = lambda==0.0
            unpinned && !vbc_final_model(model,lattice) && error("preparation removal changed the final NN model")
            row = Dict{String,Any}("branch"=>branch,"batch"=>stage_index,"stage"=>stage["name"],
                "lambda"=>lambda,"cumulative_sweeps"=>2stage_index,
                "unpinned_sweeps"=>branch=="random" ? 2stage_index : max(0,2stage_index-4),
                "status"=>"running","active_phase"=>"DMRG","final_nn_model_verified"=>unpinned,
                "model_configuration"=>KagomeDMRG._checkpoint_configuration(model,:seam,nothing;Q=6))
            push!(record["batches"],row); persist()
            function progress(event)
                if event.kind in (:sweep,:phase)
                    row["last_event"] = Dict(string(k)=>(v isa Symbol ? string(v) : v) for (k,v) in pairs(event))
                    record["worker_elapsed_seconds"] = elapsed(); persist()
                    if event.kind==:sweep
                        println(branch," ",stage["name"]," sweep=",2(stage_index-1)+event.sweep,
                            " E=",event.energy," chi=",event.maxlinkdim," elapsed=",round(elapsed();digits=1))
                        flush(stdout)
                    end
                end
            end
            row["dmrg_seconds"] = @elapsed result = run_dmrg(model,0.0;Q=6,sites,
                psi0=states[branch],seed=cfg["seed"],nsweeps=2,maxdim=stage["maxdim"],
                cutoff=0.0,noise=0.0,eigsolve_tol=1e-11,eigsolve_krylovdim=40,
                eigsolve_maxiter=20,measure_variance=false,progress_callback=progress)
            all(iszero,result.hz) || error("unexpected field")
            merge!(row,Dict("energy"=>result.energy,"energy_per_site"=>result.energy/54,
                "sz_profile"=>result.sz,"sweep_energies"=>result.sweep_energies,
                "measured_truncation_errors"=>result.max_truncation_errors,
                "settings"=>asdict(result.settings),"maxlinkdim"=>maxlinkdim(result.psi),
                "last_optimizer_energy_error"=>abs(result.energy-result.local_energy),
                "last_sweep_energy_change"=>abs(diff(result.sweep_energies)[end]),
                "column_sz"=>[sum(result.sz[9x+1:9(x+1)]) for x in 0:5],"active_phase"=>"checkpoint"))
            persist()
            row["checkpoint_save_seconds"] = @elapsed snapshot = save_checkpoint(joinpath(output,branch),result;
                baseline=result,theta_path=[0.0],status=:trial)
            row["checkpoint"] = relpath(snapshot,ROOT)
            row["checkpoint_sha256"] = Dict(p=>file_hash(joinpath(snapshot,p)) for p in ("metadata.toml","state.jls","checksums.toml"))
            persist()
            row["checkpoint_reload_seconds"] = @elapsed saved = load_checkpoint(snapshot,model;Q=6,status=:trial,
                expected_theta=0.0,expected_settings=result.settings)
            row["checkpoint_overlap_error"] = abs(abs(inner(saved.psi,result.psi))-1)
            row["checkpoint_verified"] = row["checkpoint_overlap_error"]<=1e-12 &&
                saved.diagnostics.energy==result.energy && saved.diagnostics.sz==result.sz
            row["checkpoint_verified"] || error("checkpoint reload differs")
            row["status"] = "checkpoint_saved_diagnostics_deferred"; persist()
            if stage["measure"]
                unpinned || error("structural comparison requires preparation removal")
                measure_state(result,lattice,record,row,persist)
                row["integrity_passed"] || error("observable integrity failed")
                haskey(previous,branch) && vbc_compare_precision!(row,previous[branch],54)
                previous[branch] = row
            end
            states[branch] = result.psi
            persist()
        end
        record["status"] = "completed_preparation_comparison_accuracy_separate"
    catch exception
        record["status"] = "failed_or_interrupted"
        record["exception_type"] = string(typeof(exception))
        rethrow()
    finally
        record["worker_elapsed_seconds"] = elapsed()
        record["sources_unchanged"] = all(file_hash(joinpath(ROOT,p))==h for (p,h) in sources) && file_hash(config_path)==config_hash
        record["backend_unchanged"] = KagomeDMRG._checkpoint_provenance()["source_sha256"]==record["code"]["source_sha256"]
        record["sources_unchanged"] && record["backend_unchanged"] || (record["status"]="source_changed")
        persist()
    end
    return record
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==2 || error("expected OUTPUT CONFIG.toml")
    record = vbc_main(abspath(ARGS[1]),abspath(ARGS[2]))
    println(record["status"])
    record["status"]=="completed_preparation_comparison_accuracy_separate" || exit(1)
end
