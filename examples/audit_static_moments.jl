#!/usr/bin/env julia
# Supplement the immutable chi256 records with a fresh complex H†H contraction.
const MOMENT_START = time_ns()
using KagomeDMRG, ITensors, ITensorMPS, LinearAlgebra, SHA, TOML, Dates
BLAS.set_num_threads(1)
ITensors.disable_threaded_blocksparse()

const MOMENT_ROOT = normpath(joinpath(@__DIR__, ".."))
const MOMENT_BACKEND = "outputs/source-snapshots/n27-refine-before-chirality-20260919"
const MOMENT_PARENT = "outputs/p4-campaign-nn27-q3-chi256/validation.toml"
const MOMENT_PARENT_SHA = "c68b17de275f50c6850f23b0dd3243a51f56326abb07df4bd42600b94d4bbbe1"
const MOMENT_SOURCES = ("examples/audit_static_moments.jl", "examples/run_static_moments.py",
    "examples/run_static27.py", "examples/run_static18.py")
moment_elapsed() = Float64(time_ns()-MOMENT_START)/1e9
moment_path(path) = joinpath(MOMENT_ROOT,path)
moment_sha(path) = bytes2hex(sha256(read(path)))
moment_unchanged(hashes) = all(isfile(moment_path(p)) && moment_sha(moment_path(p))==h for (p,h) in hashes)

function audit_static_moments(output)
    isdir(output) && readdir(output)==["execution.toml"] || error("use fresh launcher output")
    execution = TOML.parsefile(joinpath(output,"execution.toml"))
    execution["status"]=="running" && 0<execution["wall_limit_seconds"]<=120 || error("invalid supervisor")
    Threads.nthreads()==BLAS.get_num_threads()==1 || error("one Julia/BLAS thread required")
    realpath(pkgdir(KagomeDMRG))==realpath(moment_path(MOMENT_BACKEND)) || error("archived backend required")
    sources = Dict(p=>moment_sha(moment_path(p)) for p in MOMENT_SOURCES)
    inputs = Dict(MOMENT_PARENT=>MOMENT_PARENT_SHA)
    code = KagomeDMRG._checkpoint_provenance()
    record = Dict{String,Any}("schema_version"=>1,"status"=>"running",
        "recorded_at_utc"=>string(now(UTC)),"N"=>27,"Q"=>3,"theta"=>0.0,
        "scope"=>"separate_complex_second_moment_audit_of_saved_trials",
        "parent_record"=>MOMENT_PARENT,"parent_record_sha256"=>MOMENT_PARENT_SHA,
        "parent_records_modified"=>false,"new_DMRG"=>false,"new_ED"=>false,
        "hamiltonian_reference"=>"fresh_MPO_from_same_validated_builder_not_independent_ED",
        "original_imaginary_moment"=>"not_saved_in_V1_not_retroactively_reconstructed_as_original_measurement",
        "variance_interpretation"=>"H_dispersion_not_ground_energy_error_bound",
        "norm_charge_overlap_tolerance"=>1e-12,
        "energy_relative_tolerance"=>1e-10,
        "backend_project"=>MOMENT_BACKEND,
        "backend_snapshot_origin"=>"48d15e2b566b4d4b50a6c5295b84348cc5f30f8d",
        "git_revision_semantics"=>"execution_directory_git_lookup_not_snapshot_origin",
        "code"=>code,"analysis_source_sha256"=>sources,"input_sha256"=>inputs,
        "runtime"=>merge(KagomeDMRG._checkpoint_runtime(),Dict("julia_threads"=>1,"blas_threads"=>1)),
        "wall_limit_seconds"=>execution["wall_limit_seconds"],"points"=>Dict{String,Any}[])
    persist() = begin
        record["worker_elapsed_seconds"] = moment_elapsed()
        KagomeDMRG._write_flux_record(joinpath(output,"validation.toml"),record)
    end
    persist()
    try
        moment_unchanged(inputs) || error("pinned parent record mismatch")
        parent = TOML.parsefile(moment_path(MOMENT_PARENT))
        parent["status"]=="completed_bounded_case_accuracy_separate" || error("parent run incomplete")
        lattice = kagome_cylinder(3,3)
        configuration = KagomeDMRG._checkpoint_configuration(lattice,:seam,nothing;Q=3)
        parent["configuration"]==configuration || error("parent model mismatch")
        parent["code"]["source_sha256"]==code["source_sha256"] || error("parent backend mismatch")
        record["configuration"] = configuration
        for p in MOMENT_SOURCES
            destination = joinpath(output,"analysis-sources",p)
            mkpath(dirname(destination))
            cp(moment_path(p),destination)
        end
        for sweeps in (6,8)
            old = only(filter(r->r["cumulative_sweeps"]==sweeps,parent["batches"]))
            old["integrity_passed"] && old["checkpoint_verified"] || error("parent trial integrity missing")
            checkpoint = old["checkpoint"]
            for name in ("metadata.toml","state.jls","checksums.toml")
                p = joinpath(checkpoint,name)
                inputs[p] = moment_sha(moment_path(p))
            end
            inputs[joinpath(checkpoint,"metadata.toml")]==old["checkpoint_metadata_sha256"] || error("metadata hash mismatch")
            inputs[joinpath(checkpoint,"state.jls")]==old["checkpoint_payload_sha256"] || error("payload hash mismatch")
            row = Dict{String,Any}("cumulative_sweeps"=>sweeps,"checkpoint"=>checkpoint,
                "status"=>"running","active_phase"=>"strict_checkpoint_load")
            push!(record["points"],row)
            persist()
            meta = TOML.parsefile(moment_path(joinpath(checkpoint,"metadata.toml")))
            settings = KagomeDMRG._checkpoint_settings(meta["settings"])
            row["load_seconds"] = @elapsed saved = load_checkpoint(moment_path(checkpoint),lattice;
                Q=3,gauge=:seam,status=:trial,expected_theta=0.0,expected_settings=settings)
            saved.settings.measure_variance==false || error("expected original variance-free solver checkpoint")
            saved.diagnostics.energy==old["energy"] || error("checkpoint and parent energy differ")
            original_psi = deepcopy(saved.psi)
            row["norm_before"] = norm(saved.psi)
            row["active_phase"] = "fresh_MPO_and_energy"
            persist()
            row["MPO_seconds"] = @elapsed H = twisted_exchange_mpo(saved.sites,lattice,0.0;gauge=:seam,hz=saved.hz)
            row["energy_seconds"] = @elapsed energy = inner(saved.psi',H,saved.psi)
            row["saved_energy"] = saved.diagnostics.energy
            row["remeasured_energy_real"] = real(energy)
            row["remeasured_energy_imag"] = imag(energy)
            row["energy_difference"] = abs(energy-saved.diagnostics.energy)
            row["active_phase"] = "complex_HdaggerH"
            persist()
            row["moment_seconds"] = @elapsed moment = inner(H,saved.psi,H,saved.psi)
            variance = real(moment)-saved.diagnostics.energy^2
            floor = 100eps(Float64)*max(1.0,saved.diagnostics.energy^2)
            row["HdaggerH_real"] = real(moment)
            row["HdaggerH_imag"] = imag(moment)
            row["variance"] = variance
            row["variance_energy_convention"] = "saved_expectation_energy_after_fresh_energy_recheck"
            row["variance_using_remeasured_energy"] = real(moment)-real(energy)^2
            row["original_variance"] = old["variance"]
            row["variance_difference"] = variance-old["variance"]
            row["variance_roundoff_scale"] = floor
            row["variance_difference_tolerance"] = max(floor,old["variance_roundoff_scale"])
            row["energy_difference_tolerance"] = 1e-10*max(1,abs(saved.diagnostics.energy))
            row["norm_after"] = norm(saved.psi)
            row["overlap_change_error"] = abs(abs(inner(original_psi,saved.psi))-1)
            row["charge_error"] = abs(sum(sz_profile(saved.psi))-1.5)
            checks = Dict("finite_moment"=>isfinite(moment),
                "moment_imaginary"=>abs(imag(moment))<=floor,
                "variance_nonnegative_within_roundoff"=>isfinite(variance) && variance>=-floor,
                "original_variance_agreement"=>abs(row["variance_difference"])<=row["variance_difference_tolerance"],
                "remeasured_energy_agreement"=>isfinite(energy) && row["energy_difference"]<=row["energy_difference_tolerance"],
                "norm_before"=>abs(row["norm_before"]-1)<=1e-12,
                "norm_after"=>abs(row["norm_after"]-1)<=1e-12,
                "state_unchanged"=>row["overlap_change_error"]<=1e-12,
                "charge"=>row["charge_error"]<=1e-12)
            row["checks"] = checks
            row["passed"] = all(values(checks))
            row["active_phase"] = "completed"
            row["status"] = row["passed"] ? "passed_saved_state_moment_audit" : "failed_moment_audit"
            persist()
            row["passed"] || error("saved-state moment audit failed")
        end
        record["status"] = "passed_saved_state_moment_audit"
    catch exception
        record["status"] = exception isa InterruptException ? "interrupted" : "failed"
        record["exception_type"] = string(nameof(typeof(exception)))
        persist()
        rethrow()
    finally
        try
            record["sources_unchanged"] = moment_unchanged(sources)
            record["inputs_unchanged"] = moment_unchanged(inputs)
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
    length(ARGS)==1 || error("usage: audit_static_moments.jl OUTPUT")
    result = audit_static_moments(abspath(ARGS[1]))
    println(result["status"])
    result["status"]=="passed_saved_state_moment_audit" || exit(1)
end
