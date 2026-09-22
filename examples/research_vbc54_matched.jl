#!/usr/bin/env julia
# Keep the proven backend and parent environment unchanged. Only the research
# orchestration differs; every child starts directly from its eight-sweep parent.
include("research_vbc54.jl")

const MATCHED_SOURCES = ("examples/research_vbc54_matched.jl",
    "examples/run_vbc54_matched.py", VBC_SOURCES...)

function matched_config(path)
    cfg = TOML.parsefile(path)
    Set(keys(cfg))==Set(["schema_version","case_id","Lx","Ly","Q","seed",
        "branches","maxdims","additional_sweeps","wall_seconds","question","parents","precision"]) || error("unknown configuration keys")
    cfg["schema_version"]==1 && (cfg["Lx"],cfg["Ly"],cfg["Q"],cfg["seed"])==(6,3,6,11) || error("expected N54 Q6 seed11")
    cfg["branches"]==["random","windmill"] && cfg["maxdims"]==[128,256] &&
        cfg["additional_sweeps"]==2 && cfg["wall_seconds"]==3600 &&
        cfg["precision"]==VBC_PRECISION || error("changed comparison protocol")
    length(cfg["parents"])==2 && [p["branch"] for p in cfg["parents"]]==cfg["branches"] || error("parent branches mismatch")
    for p in cfg["parents"]
        Set(keys(p))==Set(["branch","checkpoint","record","record_sha256","metadata_sha256"]) || error("unknown parent keys")
    end
    return cfg
end

function matched_parent(p,lattice)
    evidence_path, path = rootpath(p["record"]),rootpath(p["checkpoint"])
    file_hash(evidence_path)==p["record_sha256"] || error("parent record hash mismatch")
    evidence = TOML.parsefile(evidence_path)
    evidence["status"]=="completed_preparation_comparison_accuracy_separate" &&
        evidence["sources_unchanged"] && evidence["backend_unchanged"] || error("incomplete parent research record")
    evidence["configuration"]==KagomeDMRG._checkpoint_configuration(lattice,:seam,nothing;Q=6) || error("wrong parent Hamiltonian")
    previous = only(filter(r->r["branch"]==p["branch"] && r["stage"]=="relax",evidence["batches"]))
    rootpath(previous["checkpoint"])==path && previous["cumulative_sweeps"]==8 &&
        previous["lambda"]==0 && previous["final_nn_model_verified"] &&
        previous["integrity_passed"] && previous["checkpoint_verified"] || error("wrong parent stage or integrity")
    hashes = Dict(f=>file_hash(joinpath(path,f)) for f in ("metadata.toml","state.jls","checksums.toml"))
    hashes==previous["checkpoint_sha256"] && hashes["metadata.toml"]==p["metadata_sha256"] || error("parent snapshot hash mismatch")
    cfg = Dict("parent_checkpoint"=>p["checkpoint"],"parent_record"=>p["record"],
        "parent_record_sha256"=>p["record_sha256"],"parent_metadata_sha256"=>p["metadata_sha256"],
        "initial_completed_sweeps"=>8,"Q"=>6)
    holder = Dict{String,Any}()
    saved,previous,pins = parent_state(cfg,lattice,holder)
    saved.diagnostics.energy==previous["energy"] && saved.diagnostics.sz==previous["sz_profile"] &&
        asdict(saved.settings)==previous["settings"] && saved.settings.nsweeps==2 &&
        all(==(128),saved.settings.maxdim) || error("parent diagnostics/settings mismatch")
    info = holder["parent"]
    info["checkpoint_sha256"] = hashes
    info["unpinned_sweeps"] = previous["unpinned_sweeps"]
    return saved,previous,info,pins
end

function matched_main(output,config_path)
    isdir(output) && readdir(output)==["execution.toml"] || error("use a fresh supervised output")
    cfg = matched_config(config_path)
    execution = TOML.parsefile(joinpath(output,"execution.toml"))
    execution["status"]=="running" && execution["wall_limit_seconds"]==cfg["wall_seconds"] || error("invalid supervisor")
    Threads.nthreads()==BLAS.get_num_threads()==1 || error("one Julia/BLAS thread required")
    lattice = kagome_cylinder(6,3)
    sources = Dict(p=>file_hash(joinpath(ROOT,p)) for p in MATCHED_SOURCES)
    config_hash = file_hash(config_path)
    record = Dict{String,Any}("schema_version"=>1,"case_id"=>cfg["case_id"],
        "status"=>"running","recorded_at_utc"=>string(now(UTC)),"N"=>54,"Q"=>6,"theta"=>0.0,
        "config"=>cfg,"config_sha256"=>config_hash,"analysis_source_sha256"=>sources,
        "configuration"=>KagomeDMRG._checkpoint_configuration(lattice,:seam,nothing;Q=6),
        "code"=>KagomeDMRG._checkpoint_provenance(),
        "runtime"=>merge(KagomeDMRG._checkpoint_runtime(),Dict("julia_threads"=>1,"blas_threads"=>1)),
        "wall_limit_seconds"=>cfg["wall_seconds"],"integrity_limits"=>INTEGRITY_LIMITS,
        "comparison_scope"=>"finite_additional_sweeps_from_identical_parent_not_converged_chi_extrapolation",
        "measurement_frequency"=>"full_diagnostics_once_per_child_after_two_sweeps_and_checkpoint",
        "case_order"=>"chi128_random_windmill_then_chi256_random_windmill",
        "variance_interpretation"=>"H_dispersion_not_ground_energy_error_bound",
        "parents"=>Dict{String,Any}(),"batches"=>Dict{String,Any}[])
    persist() = KagomeDMRG._write_flux_record(joinpath(output,"validation.toml"),record)
    persist()
    for p in MATCHED_SOURCES
        destination = joinpath(output,"analysis-sources",p)
        mkpath(dirname(destination)); cp(joinpath(ROOT,p),destination)
    end
    cp(config_path,joinpath(output,"analysis-sources/config.toml"))
    parents,previous = Dict{String,Any}(),Dict{String,Any}()
    parent_pins = Dict{String,String}()
    try
        for p in cfg["parents"]
            branch = p["branch"]
            parents[branch],previous[branch],record["parents"][branch],pins = matched_parent(p,lattice)
            merge!(parent_pins,pins)
            persist()
        end
        parents["random"].sites==parents["windmill"].sites || error("parents use different physical site indices")
        for chi in cfg["maxdims"], branch in cfg["branches"]
            parent = parents[branch]
            row = Dict{String,Any}("branch"=>branch,"chi"=>chi,"cumulative_sweeps"=>10,
                "unpinned_sweeps"=>previous[branch]["unpinned_sweeps"]+2,
                "parent_checkpoint"=>record["parents"][branch]["checkpoint"],
                "parent_record_sha256"=>record["parents"][branch]["record_sha256"],
                "status"=>"running","active_phase"=>"DMRG","lambda"=>0.0,"final_nn_model_verified"=>true,
                "model_configuration"=>record["configuration"])
            push!(record["batches"],row); persist()
            function progress(event)
                if event.kind in (:sweep,:phase)
                    row["last_event"] = Dict(string(k)=>(v isa Symbol ? string(v) : v) for (k,v) in pairs(event))
                    record["worker_elapsed_seconds"] = elapsed(); persist()
                    if event.kind==:sweep
                        println(branch," chi=",chi," sweep=",8+event.sweep," E=",event.energy,
                            " retained=",event.maxlinkdim," elapsed=",round(elapsed();digits=1))
                        flush(stdout)
                    end
                end
            end
            # run_dmrg copies psi0 and rebuilds the MPO/environments. The chi256
            # child receives the original parent, never the chi128 result.
            row["dmrg_seconds"] = @elapsed result = run_dmrg(lattice,0.0;Q=6,
                sites=parent.sites,psi0=parent.psi,seed=cfg["seed"],initial_linkdim=4,
                nsweeps=2,maxdim=[chi,chi],cutoff=0.0,noise=0.0,eigsolve_tol=1e-11,
                eigsolve_krylovdim=40,eigsolve_maxiter=20,measure_variance=false,progress_callback=progress)
            all(iszero,result.hz) && flux(result.psi)==QN("Sz",6) &&
                all(i->siteind(result.psi,i)==parent.sites[i],eachindex(parent.sites)) &&
                all(i->eltype(result.psi[i])<:Complex,eachindex(result.psi)) || error("child state contract failed")
            row["parent_profile_unchanged_error"] = maximum(abs.(sz_profile(parent.psi)-previous[branch]["sz_profile"]))
            row["parent_profile_unchanged_error"]<=1e-12 || error("parent mutated by solve")
            merge!(row,Dict("energy"=>result.energy,"energy_per_site"=>result.energy/54,
                "sz_profile"=>result.sz,"sweep_energies"=>result.sweep_energies,
                "measured_truncation_errors"=>result.max_truncation_errors,"settings"=>asdict(result.settings),
                "maxlinkdim"=>maxlinkdim(result.psi),"last_optimizer_energy_error"=>abs(result.energy-result.local_energy),
                "last_sweep_energy_change"=>abs(diff(result.sweep_energies)[end]),
                "column_sz"=>[sum(result.sz[9x+1:9(x+1)]) for x in 0:5],"active_phase"=>"checkpoint"))
            persist()
            row["checkpoint_save_seconds"] = @elapsed snapshot = save_checkpoint(joinpath(output,branch*"-chi"*string(chi)),
                result;baseline=result,theta_path=[0.0],status=:trial)
            row["checkpoint"] = relpath(snapshot,ROOT)
            row["checkpoint_sha256"] = Dict(p=>file_hash(joinpath(snapshot,p)) for p in ("metadata.toml","state.jls","checksums.toml"))
            persist()
            row["checkpoint_reload_seconds"] = @elapsed saved = load_checkpoint(snapshot,lattice;Q=6,status=:trial,
                sites=parent.sites,expected_theta=0.0,expected_settings=result.settings)
            row["checkpoint_overlap_error"] = abs(abs(inner(saved.psi,result.psi))-1)
            row["checkpoint_verified"] = row["checkpoint_overlap_error"]<=1e-12 &&
                saved.diagnostics.energy==result.energy && saved.diagnostics.sz==result.sz
            row["checkpoint_verified"] || error("checkpoint reload differs")
            row["status"] = "checkpoint_saved_diagnostics_deferred"; persist()
            measure_state(result,lattice,record,row,persist)
            row["integrity_passed"] || error("observable integrity failed")
            vbc_compare_precision!(row,previous[branch],54)
            row["same_chi_comparison"] = chi==128
            row["comparison_kind"] = chi==128 ? "same_chi_same_unpinned_H_additional_two_sweeps" :
                "increased_chi_same_unpinned_H_additional_two_sweeps_from_original_parent"
            persist()
            GC.gc()
        end
        record["status"] = "completed_matched_parent_comparison_accuracy_separate"
    catch exception
        record["status"] = "failed_or_interrupted"
        record["exception_type"] = string(typeof(exception))
        rethrow()
    finally
        record["worker_elapsed_seconds"] = elapsed()
        record["sources_unchanged"] = try
            all(file_hash(joinpath(ROOT,p))==h for (p,h) in sources) && file_hash(config_path)==config_hash
        catch
            false
        end
        record["backend_unchanged"] = try
            KagomeDMRG._checkpoint_provenance()["source_sha256"]==record["code"]["source_sha256"]
        catch
            false
        end
        record["parents_unchanged"] = try
            all(file_hash(rootpath(p))==h for (p,h) in parent_pins)
        catch
            false
        end
        record["sources_unchanged"] && record["backend_unchanged"] && record["parents_unchanged"] || (record["status"]="source_or_parent_changed")
        persist()
    end
    return record
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==2 || error("expected OUTPUT CONFIG.toml")
    record = matched_main(abspath(ARGS[1]),abspath(ARGS[2]))
    println(record["status"])
    record["status"]=="completed_matched_parent_comparison_accuracy_separate" || exit(1)
end
