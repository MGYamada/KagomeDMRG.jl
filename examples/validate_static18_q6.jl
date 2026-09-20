#!/usr/bin/env julia
# One new NN sector; previously measured sectors retain their original source identity.
const WORKER_START = time_ns()
module Static18Tools
include("validate_static18.jl")
end
using KagomeDMRG, ITensors, ITensorMPS, LinearAlgebra, SHA, TOML, Dates
const ROOT = normpath(joinpath(@__DIR__, ".."))
const Q = 6
const PRIOR_PATH = "docs/research/data/p4_static18_validation.toml"
const PRIOR_SHA256 = "0d0e2c5b9ab23944babea2cd605ff5d712b23a0cb5988d5a544188bfae7ac702"
const EXTRA_SOURCES = ("examples/validate_static18_q6.jl", "examples/run_static18_q6.py",
    "examples/validate_static18.jl", "examples/run_static18.py", "examples/run_static27.py",
    "test/reference_ed.jl", "test/reference_eigensolve.jl", "test/itensor_helpers.jl")
file_hash(path) = bytes2hex(sha256(read(joinpath(ROOT,path))))
extra_hashes() = Dict(path=>file_hash(path) for path in EXTRA_SOURCES)
elapsed() = (time_ns()-WORKER_START)/1e9
asdict(value) = Dict(string(k)=>v for (k,v) in pairs(value))

function load_prior(lattice)
    file_hash(PRIOR_PATH)==PRIOR_SHA256 || error("prior evidence hash mismatch")
    prior = TOML.parsefile(joinpath(ROOT,PRIOR_PATH))
    prior["status"]=="passed_finite_cluster_comparison" && prior["extra_sources_unchanged"] ||
        error("prior study did not pass")
    prior["N"]==18 && prior["target_Q"]==2 && prior["planned_Q"]==[0,2,4] ||
        error("wrong prior study")
    length(prior["runs"])==3 && sort(getindex.(prior["runs"],"Q"))==[0,2,4] ||
        error("prior charges are incomplete or duplicated")
    for path in ("test/reference_ed.jl","test/reference_eigensolve.jl","test/itensor_helpers.jl",
                 "examples/validate_static18.jl")
        prior["extra_source_sha256"][path]==file_hash(path) || error("shared reference/helper changed")
    end
    runtime = merge(KagomeDMRG._checkpoint_runtime(),Dict(
        "julia_threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads()))
    prior["runtime"]==runtime || error("prior runtime differs; explicit compatibility audit required")
    for run in prior["runs"]
        run["status"]=="passed" && run["reference_valid"] &&
            run["reference"]["converged"] && run["reference"]["boundary_level_resolved"] ||
            error("unresolved prior sector")
        cfg = KagomeDMRG._checkpoint_configuration(lattice,:seam,nothing;Q=run["Q"])
        run["configuration"]==cfg || error("prior model/geometry/charge configuration differs")
        row = last(run["attempts"])
        row["passed"] && row["checkpoint_passed"] || error("invalid prior final state")
        all(isfinite(row[k]) && row[k]<=limit for (k,limit) in prior["limits"]) ||
            error("prior numerical thresholds do not pass")
        all(isfinite,run["reference"]["values"]) || error("nonfinite prior ED energy")
    end
    return prior
end

function field_envelope(energies; target=2)
    charges = sort(collect(keys(energies)))
    haskey(energies,target) && all(isfinite,values(energies)) ||
        throw(ArgumentError("finite energies including the target are required"))
    lower = [(Q=q,h=2*(energies[target]-energies[q])/(target-q)) for q in charges if q<target]
    upper = [(Q=q,h=2*(energies[q]-energies[target])/(q-target)) for q in charges if q>target]
    !isempty(lower) && !isempty(upper) || throw(ArgumentError("both sides of the target are required"))
    lo, hi = lower[argmax(getproperty.(lower,:h))], upper[argmin(getproperty.(upper,:h))]
    return Dict("compared_Q"=>charges,"target_Q"=>target,"h_lower"=>lo.h,"h_upper"=>hi.h,
        "lower_competitor_Q"=>lo.Q,"upper_competitor_Q"=>hi.Q,"width"=>hi.h-lo.h,
        "positive_width"=>hi.h>lo.h,"lower_candidates"=>asdict.(lower),"upper_candidates"=>asdict.(upper))
end

function summarize(prior,ed,row)
    previous = Dict(run["Q"]=>run for run in prior["runs"])
    ed_energies = Dict(q=>first(run["reference"]["values"]) for (q,run) in previous)
    dmrg_energies = Dict(q=>last(run["attempts"])["energy"] for (q,run) in previous)
    old_ed,old_dmrg = field_envelope(ed_energies),field_envelope(dmrg_energies)
    ed_energies[Q],dmrg_energies[Q] = first(ed.values),row["energy"]
    new_ed,new_dmrg = field_envelope(ed_energies),field_envelope(dmrg_energies)
    h26_ed = (ed_energies[6]-ed_energies[2])/2
    h26_dmrg = (dmrg_energies[6]-dmrg_energies[2])/2
    margin = h26_ed-old_ed["h_upper"]
    # This tolerance only labels a numerical tie; it does not change any energy
    # or enlarge the interval, nor is it a rigorous spectral error bound.
    decision = margin>1e-10 ? "Q6_does_not_preempt_Q2_interval" :
        margin < -1e-10 ? (new_ed["positive_width"] ? "Q6_narrows_Q2_interval" :
            "Q6_removes_positive_width_Q2_interval") : "unresolved_near_tie"
    changes = [Dict("from_Q"=>q,"to_Q"=>6,"delta_total_sz"=>(6-q)/2,
        "site_delta_sz"=>row["sz_profile"]-last(previous[q]["attempts"])["sz_profile"],
        "column_delta_sz"=>row["column_sz"]-last(previous[q]["attempts"])["column_sz"],
        "interpretation"=>"selected_finite_cluster_states_no_edge_bulk_separation") for q in (2,4)]
    consistency = Dict{String,Any}("ed_su2_energy_order"=>all(ed_energies[a]<=ed_energies[b]+1e-10
            for (a,b) in ((0,2),(2,4),(4,6))),
        "dmrg_su2_energy_order"=>all(dmrg_energies[a]<=dmrg_energies[b]+1e-10
            for (a,b) in ((0,2),(2,4),(4,6))),
        "max_delta_total_sz_error"=>maximum(abs(sum(r["column_delta_sz"])-r["delta_total_sz"]) for r in changes))
    consistency["passed"] = consistency["ed_su2_energy_order"] && consistency["dmrg_su2_energy_order"] &&
        consistency["max_delta_total_sz_error"]<=1e-10
    return Dict("ed"=>new_ed,"dmrg"=>new_dmrg,"prior_ed"=>old_ed,"prior_dmrg"=>old_dmrg,
        "ed_energies"=>[Dict("Q"=>q,"energy"=>ed_energies[q]) for q in (0,2,4,6)],
        "dmrg_energies"=>[Dict("Q"=>q,"energy"=>dmrg_energies[q]) for q in (0,2,4,6)],
        "h26_ed"=>h26_ed,"h26_dmrg"=>h26_dmrg,"h26_method_difference"=>abs(h26_ed-h26_dmrg),
        "h26_minus_previous_upper"=>margin,"decision_tolerance"=>1e-10,"decision"=>decision,
        "Q6_minus_Q2_field_energy_at_previous_upper"=>2margin,
        "h_lower_method_difference"=>abs(new_ed["h_lower"]-new_dmrg["h_lower"]),
        "h_upper_method_difference"=>abs(new_ed["h_upper"]-new_dmrg["h_upper"]),
        "magnetization_changes"=>changes,"consistency"=>consistency,
        "unexplored_nonnegative_Q"=>collect(8:2:18),
        "negative_Q"=>"not_computed_spin_reversal_mirrors_at_theta0_hz0_for_nonnegative_field",
        "bulk_plateau"=>"not_established_two_columns_no_interior",
        "method_discrepancy"=>"observed_agreement_not_a_rigorous_error_bound")
end

function main(output)
    isdir(output) && readdir(output)==["execution.toml"] || error("use the launcher with a fresh output")
    execution = TOML.parsefile(joinpath(output,"execution.toml"))
    execution["status"]=="running" && 0<execution["wall_limit_seconds"]<=300 || error("invalid supervisor")
    Threads.nthreads()==BLAS.get_num_threads()==1 || error("one Julia/BLAS thread required")
    lattice = kagome_cylinder(2,3)
    prior = load_prior(lattice)
    before = extra_hashes()
    code = KagomeDMRG._checkpoint_provenance()
    old_sources,current_sources = prior["code"]["source_sha256"],code["source_sha256"]
    changes = sort([p for p in union(keys(old_sources),keys(current_sources))
        if get(old_sources,p,nothing)!=get(current_sources,p,nothing)])
    run = Dict{String,Any}("Q"=>Q,"theta"=>0.0,"status"=>"running",
        "configuration"=>KagomeDMRG._checkpoint_configuration(lattice,:seam,nothing;Q),
        "attempts"=>Dict{String,Any}[])
    record = Dict{String,Any}("schema_version"=>1,"status"=>"running",
        "recorded_at_utc"=>string(now(UTC)),"scope"=>"NN_N18_Q6_skipped_sector_check",
        "model"=>"nearest_neighbor_isotropic_Heisenberg_J1=1_J2=J3=hz=theta=0",
        "N"=>18,"target_Q"=>2,"new_Q"=>[Q],"compared_Q"=>[0,2,4,6],
        "code"=>code,"extra_source_sha256"=>before,"runtime"=>prior["runtime"],
        "prior_evidence"=>Dict("path"=>PRIOR_PATH,"sha256"=>PRIOR_SHA256,
            "recorded_at_utc"=>prior["recorded_at_utc"],"code"=>prior["code"],
            "runtime"=>prior["runtime"],"solver"=>prior["solver"],"ed_solver"=>prior["ed_solver"],
            "original_Q"=>[0,2,4],"configuration_matches"=>true,"shared_helpers_match"=>true,
            "different_production_source_paths"=>changes,
            "prior_energies"=>"quoted_original_measurements_not_recomputed_or_relabeled",
            "prior_checkpoints"=>"not_loaded_or_modified"),
        "solver"=>asdict(Static18Tools.SOLVER),"ed_solver"=>asdict(Static18Tools.ED_SOLVER),
        "limits"=>Static18Tools.LIMITS,"ground_cluster_tolerance"=>Static18Tools.CLUSTER_TOLERANCE,
        "retry_policy"=>"one_additional_six_sweeps_only_if_ED_comparison_fails_within_wall_limit",
        "wall_limit_seconds"=>execution["wall_limit_seconds"],"execution_record"=>"execution.toml",
        "time_limit_enforcement"=>"external_including_startup_and_single_solves",
        "phase_identification"=>"not_attempted","quantized_pump"=>"not_measured","runs"=>[run])
    persist() = KagomeDMRG._write_flux_record(joinpath(output,"validation.toml"),record)
    function phase!(name)
        run["active_phase"] = name
        record["worker_elapsed_seconds"] = elapsed()
        persist()
        println("Q=6 phase=",name," elapsed=",round(elapsed();digits=2))
        flush(stdout)
    end
    persist()
    try
        phase!("independent_ED")
        ed_seconds = @elapsed ed = Static18Tools.reference_eigensystem(2,3,0.0;Q,Static18Tools.ED_SOLVER...)
        cluster = findall(e->e-first(ed.values)<Static18Tools.CLUSTER_TOLERANCE,ed.values)
        resolved = !isempty(cluster) && length(cluster)<length(ed.values)
        run["reference"] = Dict("values"=>ed.values,"residual_norms"=>ed.residual_norms,
            "orthogonality_error"=>ed.orthogonality_error,"converged"=>ed.converged,
            "basis_dimension"=>length(ed.basis),"settings"=>asdict(ed.settings),
            "solver_info"=>asdict(ed.solver_info),"ground_cluster_size"=>length(cluster),
            "boundary_level_resolved"=>resolved,"ground_space_completeness"=>"partial_spectrum_not_proven",
            "elapsed_seconds"=>ed_seconds,"payload"=>Static18Tools.save_reference(output,Q,ed))
        resolved && (run["reference"]["gap_above_ground_cluster"] = ed.values[length(cluster)+1]-first(ed.values))
        run["reference_valid"] = ed.converged && resolved
        persist()
        if !run["reference_valid"]
            run["status"] = record["status"] = "unresolved_reference"
            return record
        end
        previous = nothing
        for attempt in 1:2
            phase!(attempt==1 ? "DMRG_6_sweeps" : "DMRG_additional_6_sweeps")
            function progress(event)
                if event.kind==:sweep
                    run["last_sweep"] = asdict(event)
                    run["last_sweep"]["kind"] = string(event.kind)
                    run["last_sweep"]["attempt"] = attempt
                    persist()
                    println("Q=6 attempt=",attempt," sweep=",event.sweep," energy=",event.energy,
                        " maxdim=",event.maxlinkdim)
                    flush(stdout)
                end
            end
            dmrg_seconds = @elapsed result = previous===nothing ?
                run_dmrg(lattice,0.0;Q,progress_callback=progress,Static18Tools.SOLVER...) :
                run_dmrg(lattice,0.0;Q,sites=previous.sites,psi0=previous.psi,
                    progress_callback=progress,Static18Tools.SOLVER...)
            phase!("checkpoint")
            snapshot = save_checkpoint(output,result;baseline=result,theta_path=[0.0],status=:trial)
            run["latest_checkpoint"] = relpath(snapshot,output)
            persist()
            phase!("diagnostics")
            diagnostics_seconds = @elapsed row = Static18Tools.measure_state(result,ed)
            row["attempt"] = attempt
            row["total_sweeps"] = 6attempt
            row["dmrg_seconds"] = dmrg_seconds
            row["diagnostics_seconds"] = diagnostics_seconds
            row["checkpoint"] = relpath(snapshot,output)
            push!(run["attempts"],row)
            phase!("reload_verification")
            loaded = load_checkpoint(snapshot,lattice;Q,status=:trial,expected_theta=0.0)
            row["checkpoint_overlap_error"] = abs(abs(inner(loaded.psi,result.psi))-1)
            row["checkpoint_passed"] = loaded.Q==Q && row["checkpoint_overlap_error"]<1e-12
            row["checkpoint_passed"] || error("checkpoint reload failed")
            persist()
            println("Q=6 attempt=",attempt," residual=",row["residual_norm"]," passed=",row["passed"])
            flush(stdout)
            row["passed"] && break
            previous = result
        end
        row = last(run["attempts"])
        run["status"] = row["passed"] ? "passed" : "unresolved_dmrg_accuracy"
        if row["passed"]
            record["sector_comparison"] = summarize(prior,ed,row)
            record["status"] = record["sector_comparison"]["consistency"]["passed"] ?
                "passed_finite_cluster_comparison" : "unresolved_sector_consistency"
        else
            record["status"] = run["status"]
        end
        phase!("finished")
    catch exception
        run["status"] = record["status"] = exception isa InterruptException ? "interrupted" : "failed"
        record["exception_type"] = string(nameof(typeof(exception)))
        rethrow()
    finally
        record["worker_elapsed_seconds"] = elapsed()
        record["extra_sources_unchanged"] = extra_hashes()==before
        record["prior_evidence_unchanged"] = file_hash(PRIOR_PATH)==PRIOR_SHA256
        record["extra_sources_unchanged"] && record["prior_evidence_unchanged"] ||
            (record["status"]="source_or_evidence_changed")
        persist()
    end
    return record
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==1 || error("usage: run_static18_q6.py fresh-output-directory")
    record = main(abspath(only(ARGS)))
    record["status"]=="passed_finite_cluster_comparison" || error("Q6 study unresolved; inspect validation.toml")
end
