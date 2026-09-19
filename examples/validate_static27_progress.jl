#!/usr/bin/env julia
# One new N27 point; save a completed trial before separate diagnostics.
const WORKER_START = time_ns()
using KagomeDMRG, ITensors, ITensorMPS, LinearAlgebra, SHA, TOML, Dates
BLAS.set_num_threads(1)
ITensors.disable_threaded_blocksparse()
const ROOT = normpath(joinpath(@__DIR__, ".."))
const SOLVER = (; seed=11, nsweeps=2, maxdim=256, cutoff=0.0, noise=0.0,
    eigsolve_tol=1e-11, eigsolve_krylovdim=40, eigsolve_maxiter=20,
    measure_variance=false)
const N = 27
const Q = 3
const EXTRA_SOURCES = ("examples/validate_static27_progress.jl",
    "examples/run_static27_progress.py", "examples/run_static27.py", "examples/run_static18.py")
asdict(value) = Dict(string(k)=>v for (k,v) in pairs(value))
elapsed() = (time_ns()-WORKER_START)/1e9
extra_hashes() = Dict(p=>bytes2hex(sha256(read(joinpath(ROOT,p)))) for p in EXTRA_SOURCES)
column_sz(profile) = [sum(profile[(9x+1):(9x+9)]) for x in 0:2]
real_rows(matrix) = [real.(collect(matrix[i,:])) for i in axes(matrix,1)]
imag_rows(matrix) = [imag.(collect(matrix[i,:])) for i in axes(matrix,1)]

function main(output)
    isdir(output) && readdir(output)==["execution.toml"] || error("use launcher with fresh output")
    execution = TOML.parsefile(joinpath(output,"execution.toml"))
    execution["status"] == "running" || error("launcher is not running")
    Threads.nthreads() == BLAS.get_num_threads() == 1 || error("one Julia/BLAS thread required")
    lattice = kagome_cylinder(3,3)
    before = extra_hashes()
    record = Dict{String,Any}("schema_version"=>1,"status"=>"running",
        "recorded_at_utc"=>string(now(UTC)),"scope"=>"N27_two_sweep_progress_and_separate_diagnostics",
        "N"=>N,"Q"=>Q,"theta"=>0.0,"solver"=>asdict(SOLVER),
        "configuration"=>KagomeDMRG._checkpoint_configuration(lattice,:seam,nothing;Q),
        "code"=>KagomeDMRG._checkpoint_provenance(),"extra_source_sha256"=>before,
        "runtime"=>merge(KagomeDMRG._checkpoint_runtime(),Dict("julia_threads"=>1,"blas_threads"=>1)),
        "wall_limit_seconds"=>execution["wall_limit_seconds"],
        "time_limit_enforcement"=>"external_launcher_including_startup_all_phases_and_termination_grace",
        "point_budget"=>1,"automatic_retry"=>false,
        "accuracy_status"=>"unestablished_two_sweeps_single_chi_single_seed",
        "independent_ED"=>"not_attempted_reference_limited_to_N18",
        "local_Krylov_convergence_flags"=>"not_exposed_by_backend",
        "field_interval"=>"not_measured_only_Q3","bulk_plateau"=>"not_established",
        "phase_identification"=>"not_attempted","axial_pump"=>"not_measured_theta_zero_only",
        "central_sites"=>collect(10:18),"cuts"=>[9,18],
        "central_region_status"=>"one_interior_column_not_established_bulk",
        "phases"=>Dict{String,Any}[])
    persist() = KagomeDMRG._write_flux_record(joinpath(output,"validation.toml"),record)
    progress = Dict{String,Any}("schema_version"=>1,"status"=>"running",
        "status_semantics"=>"last_atomic_snapshot_consult_execution_for_process_exit",
        "energy_semantics"=>"bond_and_sweep_energies_are_local_optimizer_values_not_final_MPO_expectations",
        "events"=>Dict{String,Any}[])
    callback_seconds = Ref(0.0)
    function progress_callback(event)
        started = time_ns()
        item = Dict{String,Any}(string(k)=>(v isa Symbol ? string(v) : v) for (k,v) in pairs(event))
        item["worker_elapsed_seconds"] = elapsed()
        push!(progress["events"],item)
        progress["last_event"] = item
        progress["callback_seconds_before_this_write"] = callback_seconds[]
        KagomeDMRG._write_flux_record(joinpath(output,"progress.toml"),progress)
        if event.kind == :sweep
            println("N27 completed sweep=",event.sweep," local_energy=",event.energy,
                " truncation=",event.max_truncation_error," elapsed=",round(elapsed();digits=2))
            flush(stdout)
        end
        callback_seconds[] += (time_ns()-started)/1e9
        return nothing
    end
    function phase(f,name)
        record["active_phase"] = name
        item = Dict{String,Any}("name"=>name,"status"=>"running",
            "started_worker_seconds"=>elapsed(),"start_maxrss_bytes"=>Sys.maxrss())
        push!(record["phases"],item)
        persist()
        println("N27 phase=",name," elapsed=",round(elapsed();digits=2))
        flush(stdout)
        value = f()
        item["finished_worker_seconds"] = elapsed()
        item["wall_seconds"] = item["finished_worker_seconds"]-item["started_worker_seconds"]
        item["end_maxrss_bytes"] = Sys.maxrss()
        item["status"] = "completed"
        persist()
        return value
    end
    persist()
    try
        result = phase("two_sweeps_without_variance") do
            run_dmrg(lattice,0.0;Q,SOLVER...,progress_callback)
        end
        progress["status"] = "dmrg_call_completed"
        progress["callback_seconds"] = callback_seconds[]
        KagomeDMRG._write_flux_record(joinpath(output,"progress.toml"),progress)
        record["state"] = Dict{String,Any}("energy"=>result.energy,
            "energy_per_site"=>result.energy/N,"local_energy"=>result.local_energy,
            "norm"=>norm(result.psi),"total_sz"=>sum(result.sz),"sz_profile"=>result.sz,
            "column_sz"=>column_sz(result.sz),"maxlinkdim"=>maxlinkdim(result.psi),
            "link_dimensions"=>[dim(linkind(result.psi,b)) for b in 1:(N-1)],
            "sweep_energies"=>result.sweep_energies,
            "measured_truncation_errors"=>result.max_truncation_errors,
            "variance_measured_in_dmrg"=>false,"settings"=>asdict(result.settings))
        persist()
        norm_error = abs(norm(result.psi)-1)
        charge_error = abs(sum(result.sz)-Q/2)
        record["state_integrity"] = Dict("norm_error"=>norm_error,"total_sz_error"=>charge_error,
            "passed"=>norm_error<=1e-12 && charge_error<=1e-12 && isfinite(result.energy) &&
            all(isfinite,result.sz) && length(result.sweep_energies)==SOLVER.nsweeps &&
            length(result.max_truncation_errors)==SOLVER.nsweeps &&
            all(isfinite,result.sweep_energies) &&
            all(x->isfinite(x) && 0<=x<=1,result.max_truncation_errors))
        persist()
        record["state_integrity"]["passed"] || error("completed state integrity failed")
        saved = phase("checkpoint_save_and_reload") do
            snapshot = save_checkpoint(output,result;baseline=result,theta_path=[0.0],status=:trial)
            record["checkpoint"] = relpath(snapshot,output)
            record["checkpoint_payload_bytes"] = filesize(joinpath(snapshot,"state.jls"))
            persist()
            loaded = load_checkpoint(snapshot,lattice;Q,sites=result.sites,status=:trial,
                expected_theta=0.0,expected_settings=result.settings)
            error_overlap = abs(abs(inner(loaded.psi,result.psi))-1)
            record["checkpoint_overlap_error"] = error_overlap
            record["checkpoint_verified"] = error_overlap<=1e-12 && loaded.Q==Q &&
                loaded.diagnostics.variance===nothing
            record["checkpoint_verified"] || error("checkpoint reload failed")
            loaded
        end
        phase("Schmidt_cuts_9_and_18") do
            diagnostics = Dict{String,Any}[]
            for b in (9,18)
                s = schmidt_diagnostics(result.psi,b)
                density = sum(result.sz[1:b])
                item = merge(asdict(s),Dict("density_mean_left_sz"=>density,
                    "density_error"=>abs(s.mean_left_sz-density),
                    "normalization_error"=>abs(sum(s.probabilities)-1)))
                push!(diagnostics,item)
                record["schmidt"] = diagnostics
                persist()
                item["density_error"]<=1e-10 && item["normalization_error"]<=1e-12 ||
                    error("Schmidt readout inconsistent")
            end
        end
        phase("correlations_and_bond_profile") do
            c = spin_correlations(result.psi)
            energies = [b.Jz*real(c.zz[b.i,b.j])+b.Jxy*real(c.pm[b.i,b.j]) for b in lattice.bonds]
            record["correlations"] = Dict("zz_real"=>real_rows(c.zz),"zz_imag"=>imag_rows(c.zz),
                "pm_real"=>real_rows(c.pm),"pm_imag"=>imag_rows(c.pm))
            record["bond_observables"] = Dict("energy_profile"=>energies,"energy_sum"=>sum(energies))
            errors = Dict("bond_energy_sum"=>abs(sum(energies)-result.energy),
                "zz_diagonal"=>maximum(abs.(diag(c.zz).-0.25)),
                "pm_diagonal"=>maximum(abs.(diag(c.pm).-(0.5 .+ result.sz))),
                "zz_imaginary"=>maximum(abs.(imag.(c.zz))),
                "zz_symmetry"=>maximum(abs.(c.zz-transpose(c.zz))),
                "pm_hermiticity"=>maximum(abs.(c.pm-c.pm')),
                "fixed_charge_zz_row"=>maximum(abs.(vec(sum(c.zz;dims=2))-(Q/2).*result.sz)),
                "schmidt_variance"=>maximum(abs(s["variance_left_sz"]-
                    (real(sum(c.zz[1:b,1:b]))-sum(result.sz[1:b])^2))
                    for (s,b) in zip(record["schmidt"],(9,18))))
            record["correlation_integrity_errors"] = errors
            record["correlation_integrity_passed"] = all(isfinite,c.zz) && all(isfinite,c.pm) &&
                all((isfinite(v) && v <= (k in ("bond_energy_sum","schmidt_variance") ? 1e-9 : 1e-10))
                    for (k,v) in errors)
            persist()
            record["correlation_integrity_passed"] || error("correlation readout inconsistent")
        end
        phase("variance_separate_HdaggerH_contraction") do
            norm_before = norm(result.psi)
            overlap_before = abs(inner(saved.psi,result.psi))/(norm(saved.psi)*norm_before)
            started = time_ns()
            h2 = inner(result.H,result.psi,result.H,result.psi)
            contraction_seconds = (time_ns()-started)/1e9
            variance = real(h2)-result.energy^2
            floor = 100eps(Float64)*max(1,result.energy^2)
            norm_after = norm(result.psi)
            overlap_after = abs(inner(saved.psi,result.psi))/(norm(saved.psi)*norm_after)
            record["variance"] = Dict("raw_variance"=>variance,"variance_per_site"=>variance/N,
                "HdaggerH_real"=>real(h2),"HdaggerH_imag"=>imag(h2),"roundoff_scale"=>floor,
                "contraction_seconds"=>contraction_seconds,
                "norm_before"=>norm_before,"norm_after"=>norm_after,
                "overlap_with_checkpoint_before"=>overlap_before,
                "overlap_with_checkpoint_after"=>overlap_after,
                "floor_limited"=>abs(variance)<=floor,
                "energy_std_from_variance"=>sqrt(max(0.0,variance)),
                "interpretation"=>"MPO_energy_dispersion_not_ED_residual_or_ground_energy_error")
            isfinite(variance) && variance>=-floor && abs(imag(h2))<=floor || error("invalid variance")
            abs(norm_before-1)<=1e-12 && abs(norm_after-1)<=1e-12 &&
                abs(overlap_before-1)<=1e-12 && abs(overlap_after-1)<=1e-12 &&
                flux(result.psi)==QN("Sz",Q) &&
                all(siteind(result.psi,i)==saved.sites[i] for i in 1:N) ||
                error("state changed during separate diagnostics")
        end
        record["status"] = "completed_diagnostic_pilot_accuracy_unestablished"
        record["active_phase"] = "finished"
    catch exception
        record["status"] = exception isa InterruptException ? "interrupted" : "failed"
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
    println("N27 progress study: ",record["status"])
    return record["status"]=="completed_diagnostic_pilot_accuracy_unestablished" ? 0 : 2
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==1 || error("usage: run_static27_progress.py new-output-directory")
    exit(main(abspath(only(ARGS))))
end
