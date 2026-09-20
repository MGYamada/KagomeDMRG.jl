#!/usr/bin/env julia
# A new diagnostic of immutable trial bytes, not finalization of the old worker.
const CSM_START = time_ns()
using KagomeDMRG, ITensors, ITensorMPS, LinearAlgebra, SHA, TOML, Dates
BLAS.set_num_threads(1)
ITensors.disable_threaded_blocksparse()

const CSM_ROOT = normpath(joinpath(@__DIR__, ".."))
const CSM_BACKEND = "outputs/source-snapshots/csl36-pre-env-split-20260920"
const CSM_RUN = "outputs/p3-csl36-fixed-chi-8to10-20260920"
const CSM_CHECKPOINT = joinpath(CSM_RUN,"trial/checkpoint-N8JX0z")
const CSM_PINNED = Dict(
    joinpath(CSM_RUN,"validation.toml")=>"1d55e1cd1304711abca8d490ec2fadf1d4629a38035b6af8ffe9501333584102",
    joinpath(CSM_RUN,"execution.toml")=>"3a11bec1f66185609db1834b1431824c1640634e26726bfcf7841b961a9b045b",
    joinpath(CSM_RUN,"analysis-sources/config.toml")=>"b552cb8de302c3d6c70e9ab8ccffd8d669660904dffcd67da32e2df87ce3be97",
    joinpath(CSM_CHECKPOINT,"metadata.toml")=>"a7f8f870712212e0cc27683c473edbabe0e22ec3e53931f32b2ade196ebfca45",
    joinpath(CSM_CHECKPOINT,"state.jls")=>"12b6acd3895af98dcde7d99c4c63889e03e74c6526cc60918f1d2cbdcadf9cd8",
    joinpath(CSM_CHECKPOINT,"checksums.toml")=>"174684d46681bbf949c3a9604b0433fa127d4564bf6e705253e544728d07f991")
const CSM_SOURCES = ("examples/audit_csl_saved_moment.jl", "examples/run_csl_saved_moment.py",
    "examples/run_static27.py", "examples/run_static18.py")
csm_elapsed() = Float64(time_ns()-CSM_START)/1e9
csm_path(path) = joinpath(CSM_ROOT,path)
csm_sha(path) = bytes2hex(sha256(read(path)))
csm_unchanged(hashes) = all(isfile(csm_path(p)) && csm_sha(csm_path(p))==h for (p,h) in hashes)

function audit_csl_saved_moment(output)
    isdir(output) && readdir(output)==["execution.toml"] || error("use fresh launcher output")
    execution = TOML.parsefile(joinpath(output,"execution.toml"))
    execution["status"]=="running" && 0<execution["wall_limit_seconds"]<=300 || error("invalid supervisor")
    Threads.nthreads()==BLAS.get_num_threads()==1 || error("one Julia/BLAS thread required")
    realpath(pkgdir(KagomeDMRG))==realpath(csm_path(CSM_BACKEND)) || error("archived backend required")
    sources = Dict(p=>csm_sha(csm_path(p)) for p in CSM_SOURCES)
    inputs = copy(CSM_PINNED)
    code = KagomeDMRG._checkpoint_provenance()
    runtime = KagomeDMRG._checkpoint_runtime()
    row = Dict{String,Any}("status"=>"running","active_phase"=>"verify_immutable_inputs",
        "checkpoint"=>CSM_CHECKPOINT)
    record = Dict{String,Any}("schema_version"=>1,"status"=>"running",
        "recorded_at_utc"=>string(now(UTC)),"N"=>36,"Q"=>0,"theta"=>0.0,
        "scope"=>"separate_saved_CSL_trial_complex_second_moment_supplement",
        "couplings"=>Dict("J1"=>1.0,"J2"=>0.5,"J3"=>0.5),
        "parent_record"=>joinpath(CSM_RUN,"validation.toml"),
        "parent_execution"=>joinpath(CSM_RUN,"execution.toml"),
        "parent_records_modified"=>false,"new_DMRG"=>false,"new_ED"=>false,
        "original_worker_finalized"=>false,"original_full_integrity_retroactively_claimed"=>false,
        "original_variance"=>"missing_in_timed_out_worker",
        "hamiltonian_reference"=>"fresh_MPO_from_original_validated_builder_not_independent_ED",
        "variance_interpretation"=>"H_dispersion_not_ground_energy_error_bound",
        "checkpoint_solver_variance_setting_changed"=>false,
        "accuracy_status"=>"unestablished","phase_identification"=>"not_attempted",
        "quantized_pump"=>"not_established","nonzero_flux_performed"=>false,
        "norm_charge_profile_overlap_tolerance"=>1e-12,"energy_relative_tolerance"=>1e-10,
        "backend_project"=>CSM_BACKEND,
        "git_revision_semantics"=>"execution_directory_git_lookup_not_snapshot_origin",
        "code"=>code,"analysis_source_sha256"=>sources,"input_sha256"=>inputs,
        "runtime"=>merge(runtime,Dict("julia_threads"=>1,"blas_threads"=>1)),
        "wall_limit_seconds"=>execution["wall_limit_seconds"],"measurement"=>row)
    persist() = begin
        record["worker_elapsed_seconds"] = csm_elapsed()
        KagomeDMRG._write_flux_record(joinpath(output,"validation.toml"),record)
    end
    phase!(phase) = begin
        row["active_phase"] = phase
        persist()
        println("CSL saved moment phase=",phase," elapsed=",round(csm_elapsed();digits=2))
        flush(stdout)
    end
    persist()
    try
        csm_unchanged(inputs) || error("pinned input mismatch")
        old = TOML.parsefile(csm_path(joinpath(CSM_RUN,"validation.toml")))
        stopped = TOML.parsefile(csm_path(joinpath(CSM_RUN,"execution.toml")))
        meta = TOML.parsefile(csm_path(joinpath(CSM_CHECKPOINT,"metadata.toml")))
        cfg = TOML.parsefile(csm_path(joinpath(CSM_RUN,"analysis-sources/config.toml")))
        old["status"]=="running" && old["mode"]=="prepare" || error("unexpected old worker state")
        stopped["status"]=="timed_out" && stopped["worker_exit_confirmed"]===true &&
            stopped["worker_exit_code"]!=0 || error("original worker exit not confirmed")
        old["flux_observations"]==[] && old["N"]==36 && old["Q"]==0 || error("wrong original scope")
        all(!haskey(old,k) for k in ("configuration_and_source_unchanged","parent_unchanged")) ||
            error("unexpected old finalization guards")
        original = only(old["preparation_batches"])
        original["status"]=="running" && original["active_phase"]=="variance" &&
            original["sweeps_in_this_invocation"]==2 || error("not the supported stopped variance phase")
        all(!haskey(original,k) for k in ("raw_variance","variance_source","variance_roundoff_scale",
            "second_moment_imaginary_part","second_moment_imaginary_part_measured",
            "integrity_passed","integrity_failures","change_from_parent")) || error("unexpected completed diagnostics")
        original["checkpoint_overlap_measured"]===true &&
            0<=original["checkpoint_overlap_error"]<=1e-10 || error("original reload unverified")
        Set(keys(original["measurement_seconds"]))==Set(["schmidt","bond_energies","scalar_chirality"]) ||
            error("pre-variance diagnostics incomplete")
        original["snapshot"]["path"]==CSM_CHECKPOINT &&
            joinpath(CSM_RUN,original["checkpoint"])==CSM_CHECKPOINT || error("wrong checkpoint link")
        for (name,key) in (("metadata.toml","metadata_sha256"),("state.jls","state_sha256"),
                           ("checksums.toml","checksums_file_sha256"))
            inputs[joinpath(CSM_CHECKPOINT,name)]==original["snapshot"][key] || error("snapshot hash mismatch")
        end
        old["config"]==cfg && old["config_sha256"]==inputs[joinpath(CSM_RUN,"analysis-sources/config.toml")] ||
            error("original configuration mismatch")
        cfg["mode"]=="prepare" && cfg["preparation"]==Dict("batch_sweeps"=>2,"batches"=>1) ||
            error("unexpected original preparation")
        for (p,h) in old["extra_source_sha256"]
            archived = joinpath(CSM_RUN,"analysis-sources",p)
            inputs[archived] = h
        end
        csm_unchanged(inputs) || error("original archived analysis source mismatch")
        lattice = kagome_j1j2j3_cylinder(3,4;J1=1.0,J2=0.5,J3=0.5)
        configuration = KagomeDMRG._checkpoint_configuration(lattice,:seam,nothing;Q=0)
        old["configuration"]==meta["configuration"]==configuration || error("model mismatch")
        old["code"]["source_sha256"]==meta["provenance"]["source_sha256"]==code["source_sha256"] ||
            error("original backend mismatch")
        meta["runtime"]==runtime && old["runtime"]==record["runtime"] || error("original runtime mismatch")
        original["settings"]==meta["settings"] || error("saved settings mismatch")
        settings = KagomeDMRG._checkpoint_settings(meta["settings"])
        settings.nsweeps==2 && settings.maxdim==[256] && !settings.measure_variance &&
            settings.cutoff==settings.noise==0.0 || error("unexpected trial settings")
        record["configuration"] = configuration
        record["original_worker_status"] = old["status"]
        record["original_execution_status"] = stopped["status"]
        record["original_active_phase"] = original["active_phase"]
        record["original_checkpoint_reload_overlap_error"] = original["checkpoint_overlap_error"]
        for p in CSM_SOURCES
            destination = joinpath(output,"analysis-sources",p)
            mkpath(dirname(destination))
            cp(csm_path(p),destination)
        end
        phase!("strict_checkpoint_load")
        row["load_seconds"] = @elapsed saved = load_checkpoint(csm_path(CSM_CHECKPOINT),lattice;
            Q=0,gauge=:seam,status=:trial,expected_theta=0.0,expected_settings=settings)
        saved.theta_path==[0.0] && saved.diagnostics.energy==original["energy"] &&
            saved.diagnostics.sz==original["sz_profile"] || error("saved diagnostics mismatch")
        original_psi = deepcopy(saved.psi)
        row["norm_before"] = norm(saved.psi)
        phase!("fresh_MPO_and_energy")
        row["MPO_seconds"] = @elapsed H = twisted_exchange_mpo(saved.sites,lattice,0.0;gauge=:seam,hz=saved.hz)
        row["energy_seconds"] = @elapsed energy = inner(saved.psi',H,saved.psi)
        row["saved_energy"] = saved.diagnostics.energy
        row["remeasured_energy_real"] = real(energy)
        row["remeasured_energy_imag"] = imag(energy)
        row["energy_difference"] = abs(energy-saved.diagnostics.energy)
        row["energy_difference_tolerance"] = 1e-10*max(1,abs(saved.diagnostics.energy))
        phase!("complex_HdaggerH")
        row["moment_seconds"] = @elapsed moment = inner(H,saved.psi,H,saved.psi)
        variance = real(moment)-saved.diagnostics.energy^2
        floor = 100eps(Float64)*max(1.0,saved.diagnostics.energy^2)
        row["HdaggerH_real"] = real(moment)
        row["HdaggerH_imag"] = imag(moment)
        row["raw_variance"] = variance
        row["variance_energy_convention"] = "saved_expectation_energy_after_fresh_energy_recheck"
        row["variance_using_remeasured_energy"] = real(moment)-real(energy)^2
        row["variance_roundoff_scale"] = floor
        phase!("norm_charge_profile_state_recheck")
        row["norm_after"] = norm(saved.psi)
        row["overlap_change_error"] = abs(abs(inner(original_psi,saved.psi))-1)
        sz = sz_profile(saved.psi)
        row["remeasured_sz_profile"] = sz
        row["total_sz_error"] = abs(sum(sz))
        row["profile_max_difference"] = maximum(abs.(sz-original["sz_profile"]))
        # Fresh tensor differences detect changes without discarding a global phase.
        row["max_tensor_change"] = maximum(norm(saved.psi[k]-original_psi[k]) for k in eachindex(saved.psi))
        row["checks"] = Dict("finite_moment"=>isfinite(moment),
            "moment_imaginary"=>abs(imag(moment))<=floor,
            "variance_nonnegative_within_roundoff"=>isfinite(variance) && variance>=-floor,
            "remeasured_energy_agreement"=>isfinite(energy) && row["energy_difference"]<=row["energy_difference_tolerance"],
            "norm_before"=>abs(row["norm_before"]-1)<=1e-12,
            "norm_after"=>abs(row["norm_after"]-1)<=1e-12,
            "state_overlap"=>row["overlap_change_error"]<=1e-12,
            "state_tensors_unchanged"=>row["max_tensor_change"]==0.0,
            "charge"=>row["total_sz_error"]<=1e-12,
            "profile_agreement"=>all(isfinite,sz) && row["profile_max_difference"]<=1e-12)
        row["passed"] = all(values(row["checks"]))
        row["active_phase"] = "completed"
        row["status"] = row["passed"] ? "passed_saved_state_moment_audit" : "failed_moment_audit"
        record["status"] = row["status"]
    catch exception
        record["status"] = exception isa InterruptException ? "interrupted" : "failed"
        record["exception_type"] = string(nameof(typeof(exception)))
        persist()
        rethrow()
    finally
        try
            record["sources_unchanged"] = csm_unchanged(sources)
            record["inputs_unchanged"] = csm_unchanged(inputs)
            record["backend_unchanged"] = KagomeDMRG._checkpoint_provenance()["source_sha256"]==code["source_sha256"]
            all(record[k] for k in ("sources_unchanged","inputs_unchanged","backend_unchanged")) ||
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
    length(ARGS)==1 || error("usage: audit_csl_saved_moment.jl OUTPUT")
    result = audit_csl_saved_moment(abspath(ARGS[1]))
    println(result["status"])
    result["status"]=="passed_saved_state_moment_audit" || exit(1)
end
