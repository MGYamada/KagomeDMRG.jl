#!/usr/bin/env julia
# Stage B: one finite-chi N27 state and its resource/diagnostic record.
# This worker never materializes the full Hilbert space or runs ED.
const WORKER_START = time_ns()
using KagomeDMRG, ITensors, ITensorMPS, LinearAlgebra, SHA, TOML, Dates
BLAS.set_num_threads(1)
ITensors.disable_threaded_blocksparse()
const ROOT = normpath(joinpath(@__DIR__, ".."))
const SOLVER = (; seed=11, nsweeps=8, maxdim=256, cutoff=0.0, noise=0.0,
    eigsolve_tol=1e-11, eigsolve_krylovdim=40, eigsolve_maxiter=20,
    measure_variance=true)
const N = 27
const Q = 3
const CUT_BONDS = (9,18)
const EXTRA_SOURCES = ("examples/validate_static27.jl","examples/run_static27.py",
    "examples/run_static18.py")
const LIMITS = Dict("norm_error"=>1e-12,"total_sz_error"=>1e-12,
    "zz_diagonal_error"=>1e-10,"pm_diagonal_error"=>1e-10,
    "zz_imaginary_error"=>1e-10,"zz_symmetry_error"=>1e-10,
    "pm_hermiticity_error"=>1e-10,"fixed_charge_zz_row_error"=>1e-10,
    "bond_energy_sum_error"=>1e-9,"max_schmidt_density_error"=>1e-10,
    "max_schmidt_variance_error"=>1e-9,"max_schmidt_normalization_error"=>1e-12)
struct StudyTimeLimit <: Exception end
asdict(value) = Dict(string(k)=>v for (k,v) in pairs(value))
elapsed() = (time_ns()-WORKER_START)/1e9
extra_hashes() = Dict(p=>bytes2hex(sha256(read(joinpath(ROOT,p)))) for p in EXTRA_SOURCES)
real_rows(matrix) = [real.(collect(matrix[i,:])) for i in axes(matrix,1)]
imag_rows(matrix) = [imag.(collect(matrix[i,:])) for i in axes(matrix,1)]
column_sz(profile) = [sum(profile[(9x+1):(9x+9)]) for x in 0:2]

function main(output)
    isdir(output) || error("launch through run_static27.py with a new output directory")
    readdir(output) == ["execution.toml"] || error("output must contain only launcher execution.toml")
    execution = TOML.parsefile(joinpath(output,"execution.toml"))
    get(execution,"status","") == "running" || error("launcher is not running")
    Threads.nthreads() == BLAS.get_num_threads() == 1 || error("one Julia/BLAS thread required")
    compute_budget = execution["wall_limit_seconds"]-execution["termination_grace_seconds"]
    before = extra_hashes()
    lattice = kagome_cylinder(3,3)
    record = Dict{String,Any}("schema_version"=>1,"status"=>"running",
        "recorded_at_utc"=>string(now(UTC)),"scope"=>"static_N27_NN_resource_accuracy_pilot_stage_B",
        "model"=>"nearest_neighbor_isotropic_Heisenberg_J1=1_J2=J3=hz=theta=0",
        "N"=>N,"Q"=>Q,"theta"=>0.0,"solver"=>asdict(SOLVER),
        "configuration"=>KagomeDMRG._checkpoint_configuration(lattice,:seam,nothing;Q),
        "code"=>KagomeDMRG._checkpoint_provenance(),"extra_source_sha256"=>before,
        "runtime"=>merge(KagomeDMRG._checkpoint_runtime(),Dict(
            "julia_threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads())),
        "wall_limit_seconds"=>execution["wall_limit_seconds"],
        "time_limit_enforcement"=>"external_launcher_including_single_solves_and_startup",
        "execution_record"=>"execution.toml","integrity_limits"=>LIMITS,
        "point_budget"=>1,"automatic_retry"=>false,
        "accuracy_status"=>"unestablished_single_chi_single_seed",
        "accuracy_thresholds"=>"report_only_no_ED_reference_or_convergence_claim",
        "independent_ED"=>"not_attempted_reference_limited_to_N18",
        "local_Krylov_convergence_flags"=>"not_exposed_by_backend",
        "magnetization_field_interval"=>"not_measured_only_Q3",
        "bulk_plateau"=>"not_established","phase_identification"=>"not_attempted",
        "chirality"=>"not_measured","axial_pump"=>"not_measured_theta_zero_only",
        "cuts"=>collect(CUT_BONDS),"central_sites"=>collect(10:18),
        "central_region_status"=>"one_interior_column_not_established_bulk",
        "memory_method"=>"Sys.maxrss_process_high_water_bytes_including_JIT_and_runtime",
        "phases"=>Dict{String,Any}[])
    persist() = KagomeDMRG._write_flux_record(joinpath(output,"validation.toml"),record)
    function begin_phase!(phase)
        record["active_phase"] = phase
        remaining = max(0.0,compute_budget-elapsed())
        push!(record["phases"],Dict{String,Any}("name"=>phase,
            "started_worker_seconds"=>elapsed(),"remaining_seconds"=>remaining,
            "start_maxrss_bytes"=>Sys.maxrss(),"status"=>"running"))
        persist()
        remaining > 0 || throw(StudyTimeLimit())
        println("N27 Q3 phase=",phase," elapsed=",round(elapsed();digits=2))
        flush(stdout)
    end
    function end_phase!()
        phase = last(record["phases"])
        phase["finished_worker_seconds"] = elapsed()
        phase["wall_seconds"] = phase["finished_worker_seconds"]-phase["started_worker_seconds"]
        phase["end_maxrss_bytes"] = Sys.maxrss()
        phase["status"] = "completed"
        persist()
    end
    persist()
    try
        begin_phase!("DMRG_8_sweeps")
        result = run_dmrg(lattice,0.0;Q,SOLVER...)
        variance_floor = 100eps(Float64)*max(1,result.energy^2)
        record["state"] = Dict{String,Any}("energy"=>result.energy,
            "energy_per_site"=>result.energy/N,"local_energy"=>result.local_energy,
            "norm"=>norm(result.psi),"total_sz"=>sum(result.sz),
            "sz_profile"=>result.sz,"column_sz"=>column_sz(result.sz),
            "column_magnetization_ratio"=>column_sz(result.sz)./4.5,
            "central_mean_sz"=>sum(result.sz[10:18])/9,
            "central_sz_min"=>minimum(result.sz[10:18]),
            "central_sz_max"=>maximum(result.sz[10:18]),
            "raw_variance"=>result.variance,"variance_per_site"=>result.variance/N,
            "variance_roundoff_scale"=>variance_floor,
            "variance_floor_limited"=>abs(result.variance)<=variance_floor,
            "energy_std_from_variance"=>sqrt(max(0.0,result.variance)),
            "variance_interpretation"=>"MPO_energy_dispersion_not_independent_ED_residual_or_ground_energy_error",
            "maxlinkdim"=>maxlinkdim(result.psi),
            "link_dimensions"=>[dim(linkind(result.psi,b)) for b in 1:(N-1)],
            "sweep_energies"=>result.sweep_energies,
            "successive_sweep_energy_changes"=>diff(result.sweep_energies),
            "last_sweep_energy_change"=>abs(result.sweep_energies[end]-result.sweep_energies[end-1]),
            "measured_truncation_errors"=>result.max_truncation_errors,
            "last_sweep_max_truncation_error"=>last(result.max_truncation_errors),
            "settings"=>asdict(result.settings))
        end_phase!()

        # Preserve the completed state before the more expensive observables.
        begin_phase!("checkpoint_save_and_reload")
        snapshot = save_checkpoint(output,result;baseline=result,theta_path=[0.0],status=:trial)
        record["checkpoint"] = relpath(snapshot,output)
        persist()
        loaded = load_checkpoint(snapshot,lattice;Q,status=:trial,expected_theta=0.0)
        overlap_error = abs(abs(inner(loaded.psi,result.psi))-1)
        record["checkpoint_overlap_error"] = overlap_error
        record["checkpoint_verified"] = loaded.Q == Q && overlap_error < 1e-12
        record["checkpoint_verified"] || error("checkpoint verification failed")
        record["checkpoint_payload_bytes"] = filesize(joinpath(snapshot,"state.jls"))
        end_phase!()
        loaded = nothing

        begin_phase!("correlations_and_bond_profile")
        correlations = spin_correlations(result.psi)
        # Exactly theta=0, isotropic NN, no field: reuse the measured matrices.
        # Their bond sum is checked against the independently contracted MPO energy.
        energies = [b.Jz*real(correlations.zz[b.i,b.j]) +
                    b.Jxy*real(correlations.pm[b.i,b.j]) for b in lattice.bonds]
        record["correlations"] = Dict("indexing"=>"rows_and_columns_are_site_indices_1_to_27",
            "zz_real"=>real_rows(correlations.zz),"zz_imag"=>imag_rows(correlations.zz),
            "pm_real"=>real_rows(correlations.pm),"pm_imag"=>imag_rows(correlations.pm))
        record["bond_observables"] = Dict{String,Any}(
            "method"=>"theta_zero_NN_bond_energy_from_saved_zz_and_pm",
            "energy_profile"=>energies,"energy_sum"=>sum(energies),
            "bond_indices_within_columns"=>[[k for (k,b) in enumerate(lattice.bonds)
                if lattice.sites[b.i].x == x && lattice.sites[b.j].x == x] for x in 0:2])
        end_phase!()

        begin_phase!("Schmidt_cuts_9_and_18")
        schmidt_records = Dict{String,Any}[]
        for b in CUT_BONDS
            schmidt = schmidt_diagnostics(result.psi,b)
            density_mean = sum(result.sz[1:b])
            correlation_variance = real(sum(correlations.zz[1:b,1:b]))-density_mean^2
            push!(schmidt_records,merge(asdict(schmidt),Dict(
                "density_mean_left_sz"=>density_mean,
                "correlation_variance_left_sz"=>correlation_variance,
                "density_error"=>abs(schmidt.mean_left_sz-density_mean),
                "variance_error"=>abs(schmidt.variance_left_sz-correlation_variance),
                "normalization_error"=>abs(sum(schmidt.probabilities)-1))))
            record["schmidt"] = schmidt_records
            persist()
        end
        end_phase!()

        integrity = Dict{String,Any}(
            "norm_error"=>abs(norm(result.psi)-1),
            "total_sz_error"=>abs(sum(result.sz)-Q/2),
            "zz_diagonal_error"=>maximum(abs.(diag(correlations.zz).-0.25)),
            "pm_diagonal_error"=>maximum(abs.(diag(correlations.pm).-(0.5 .+ result.sz))),
            "zz_imaginary_error"=>maximum(abs.(imag.(correlations.zz))),
            "zz_symmetry_error"=>maximum(abs.(correlations.zz-transpose(correlations.zz))),
            "pm_hermiticity_error"=>maximum(abs.(correlations.pm-correlations.pm')),
            "fixed_charge_zz_row_error"=>maximum(abs.(vec(sum(correlations.zz;dims=2))-(Q/2).*result.sz)),
            "bond_energy_sum_error"=>abs(sum(energies)-result.energy),
            "max_schmidt_density_error"=>maximum(r["density_error"] for r in schmidt_records),
            "max_schmidt_variance_error"=>maximum(r["variance_error"] for r in schmidt_records),
            "max_schmidt_normalization_error"=>maximum(r["normalization_error"] for r in schmidt_records),
            "variance_physical_with_roundoff"=>isfinite(result.variance) && result.variance >= -variance_floor,
            "truncation_diagnostics_valid"=>length(result.max_truncation_errors)==SOLVER.nsweeps &&
                all(x->isfinite(x) && 0 <= x <= 1,result.max_truncation_errors),
            "sweep_diagnostics_valid"=>length(result.sweep_energies)==SOLVER.nsweeps &&
                all(isfinite,result.sweep_energies),
            "finite_observables"=>all(isfinite,result.sz) && all(isfinite,energies) &&
                all(isfinite,correlations.zz) && all(isfinite,correlations.pm))
        failed = sort([k for (k,limit) in LIMITS if !isfinite(integrity[k]) || integrity[k]>limit])
        for k in ("variance_physical_with_roundoff","truncation_diagnostics_valid",
                  "sweep_diagnostics_valid","finite_observables")
            integrity[k] || push!(failed,k)
        end
        integrity["failed_metrics"] = failed
        integrity["passed"] = isempty(failed)
        record["integrity"] = integrity
        record["status"] = isempty(failed) ? "completed_pilot_accuracy_unestablished" : "failed_observable_consistency"
        record["active_phase"] = "finished"
    catch exception
        record["status"] = exception isa StudyTimeLimit ? "time_budget_exceeded" :
            exception isa InterruptException ? "interrupted" : "failed"
        record["exception_type"] = string(nameof(typeof(exception)))
        if !isempty(record["phases"]) && last(record["phases"])["status"]=="running"
            last(record["phases"])["status"] = record["status"]
        end
        rethrow()
    finally
        record["worker_elapsed_seconds"] = elapsed()
        record["worker_peak_rss_bytes"] = Sys.maxrss()
        record["extra_sources_unchanged"] = extra_hashes()==before
        record["extra_sources_unchanged"] || (record["status"]="source_changed")
        persist()
    end
    println("N27 pilot: ",record["status"]," energy=",record["state"]["energy"],
        " variance=",record["state"]["raw_variance"],
        " truncation=",record["state"]["last_sweep_max_truncation_error"])
    return record["status"]=="completed_pilot_accuracy_unestablished" ? 0 : 2
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==1 || error("usage: run_static27.py new-output-directory")
    exit(main(abspath(only(ARGS))))
end
