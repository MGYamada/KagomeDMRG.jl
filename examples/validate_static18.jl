#!/usr/bin/env julia
# Stage A: finite-cluster NN magnetization-sector comparison, not a bulk plateau test.
# Use run_static18.py to enforce the wall limit even inside a single eigensolve.
const WORKER_START = time_ns()
using KagomeDMRG, ITensors, ITensorMPS, LinearAlgebra, SHA, TOML, Dates, Serialization
BLAS.set_num_threads(1)
ITensors.disable_threaded_blocksparse()
const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "test", "reference_ed.jl"))
include(joinpath(ROOT, "test", "reference_eigensolve.jl"))
include(joinpath(ROOT, "test", "itensor_helpers.jl"))

const CHARGES = (0, 2, 4)
const SOLVER = (; seed=11, nsweeps=6, maxdim=512, cutoff=0.0, noise=0.0,
    eigsolve_tol=1e-11, eigsolve_krylovdim=40, eigsolve_maxiter=20,
    measure_variance=true)
const ED_SOLVER = (; nev=4, blocksize=4, seed=5678, krylovdim=80,
    maxiter=200, tol=1e-11)
const CLUSTER_TOLERANCE = 1e-8
const LIMITS = Dict("energy_error_per_site"=>1e-8, "residual_norm"=>2e-6,
    "ground_space_leakage"=>2e-6, "max_sz_error"=>1e-6,
    "total_sz_error"=>1e-12, "state_norm_error"=>1e-12,
    "mpo_energy_error"=>1e-10, "bond_energy_sum_error"=>1e-10,
    "schmidt_density_error"=>1e-10)
const EXTRA_SOURCES = ("examples/validate_static18.jl", "examples/run_static18.py",
    "test/reference_ed.jl", "test/reference_eigensolve.jl", "test/itensor_helpers.jl")
struct StudyTimeLimit <: Exception end
extra_hashes() = Dict(p=>bytes2hex(sha256(read(joinpath(ROOT,p)))) for p in EXTRA_SOURCES)
elapsed() = (time_ns()-WORKER_START)/1e9
asdict(value) = Dict(string(k)=>v for (k,v) in pairs(value))
column_sz(profile) = [sum(profile[1:9]), sum(profile[10:18])]

function field_interval(energies)
    all(haskey(energies,Q) for Q in CHARGES) ||
        throw(ArgumentError("all three sector energies are required"))
    all(isfinite(energies[Q]) for Q in CHARGES) ||
        throw(ArgumentError("energies must be finite"))
    lower = energies[2]-energies[0]
    upper = energies[4]-energies[2]
    return Dict("h_lower"=>lower, "h_upper"=>upper, "width"=>upper-lower,
        "positive_width"=>upper > lower, "compared_Q"=>collect(CHARGES),
        "scope"=>"finite_N18_stability_only_within_compared_sectors")
end

function save_reference(output, Q, ed)
    path = joinpath(output,"reference_Q$(Q).jls")
    temporary, io = mktemp(output)
    try
        serialize(io,(; Q, theta=0.0, basis=ed.basis, values=ed.values, vectors=ed.vectors))
        close(io)
        mv(temporary,path)
    finally
        isopen(io) && close(io)
        isfile(temporary) && rm(temporary)
    end
    return Dict("file"=>basename(path), "sha256"=>bytes2hex(sha256(read(path))))
end

function measure_state(result, ed)
    v = sector_amplitudes(result.psi,result.sites,ed.basis)
    cluster = findall(e->e-first(ed.values) < CLUSTER_TOLERANCE,ed.values)
    ground = ed.vectors[:,cluster]
    projection = ground*(ground'*v)
    projection_norm = norm(projection)
    ed_profile = projection_norm > 0 ? reference_sz(projection,ed.basis,18) : Float64[]
    action = ed.H*v
    norm_squared = real(dot(v,v))
    rayleigh = real(dot(v,action))/norm_squared
    residual = norm(action-result.energy*v)
    energies = bond_energies(result.psi,result.lattice,0.0)
    schmidt = schmidt_diagnostics(result.psi,9)
    row = Dict{String,Any}("energy"=>result.energy, "ed_energy"=>first(ed.values),
        "energy_error_per_site"=>abs(result.energy-first(ed.values))/18,
        "residual_norm"=>residual, "ground_space_leakage"=>norm(v-projection),
        "ground_projection_norm"=>projection_norm,
        "max_sz_error"=>projection_norm > 0 ? maximum(abs.(result.sz-ed_profile)) : Inf,
        "total_sz_error"=>abs(sum(result.sz)-result.Q/2),
        "state_norm_squared"=>norm_squared,
        "state_norm_error"=>abs(norm(v)-1),
        "independent_rayleigh_energy"=>rayleigh,
        "independent_rayleigh_residual_norm"=>norm(action-rayleigh*v),
        "mpo_energy_error"=>abs(rayleigh-result.energy),
        "sz_profile"=>result.sz, "projected_ed_sz_profile"=>ed_profile,
        "column_sz"=>column_sz(result.sz),
        "projected_ed_column_sz"=>isempty(ed_profile) ? Float64[] : column_sz(ed_profile),
        "bond_energy_profile"=>energies,
        "bond_energy_sum_error"=>abs(sum(energies)-result.energy),
        "raw_variance"=>result.variance,
        "variance_residual_squared_error"=>abs(result.variance-residual^2),
        "variance_roundoff_scale"=>100eps(Float64)*max(1,result.energy^2),
        "maxlinkdim"=>maxlinkdim(result.psi),
        "sweep_energies"=>result.sweep_energies,
        "last_sweep_energy_change"=>abs(result.sweep_energies[end]-result.sweep_energies[end-1]),
        "measured_truncation_errors"=>result.max_truncation_errors,
        "schmidt_entropy"=>schmidt.entropy,
        "schmidt_mean_left_sz"=>schmidt.mean_left_sz,
        "schmidt_density_error"=>abs(schmidt.mean_left_sz-sum(result.sz[1:9])),
        "settings"=>asdict(result.settings))
    row["failed_metrics"] = sort([k for (k,limit) in LIMITS
        if !isfinite(row[k]) || row[k] > limit])
    row["passed"] = isempty(row["failed_metrics"])
    return row
end

function summarize!(record)
    runs = record["runs"]
    bycharge = Dict(r["Q"]=>r for r in runs)
    all(Q->haskey(bycharge,Q) && get(bycharge[Q],"reference_valid",false),CHARGES) || return
    ed_energies = Dict(Q=>first(bycharge[Q]["reference"]["values"]) for Q in CHARGES)
    summary = Dict{String,Any}("ed"=>field_interval(ed_energies),
        "farther_sectors"=>"not_explored_no_global_stability_claim",
        "bulk_plateau"=>"not_established_two_columns_no_bulk_interior",
        "neutral_bulk_gap"=>"not_established", "phase_identification"=>"not_attempted")
    record["sector_comparison"] = summary
    all(Q->get(bycharge[Q],"status","") == "passed",CHARGES) || return
    finals = Dict(Q=>last(bycharge[Q]["attempts"]) for Q in CHARGES)
    dmrg_energies = Dict(Q=>finals[Q]["energy"] for Q in CHARGES)
    summary["dmrg"] = field_interval(dmrg_energies)
    delta = Dict(Q=>abs(dmrg_energies[Q]-ed_energies[Q]) for Q in CHARGES)
    summary["observed_ed_dmrg_discrepancy"] = Dict(
        "h_lower"=>abs(summary["ed"]["h_lower"]-summary["dmrg"]["h_lower"]),
        "h_upper"=>abs(summary["ed"]["h_upper"]-summary["dmrg"]["h_upper"]),
        "width"=>abs(summary["ed"]["width"]-summary["dmrg"]["width"]),
        "lower_sum_of_energy_discrepancies"=>delta[0]+delta[2],
        "upper_sum_of_energy_discrepancies"=>delta[2]+delta[4],
        "width_sum_of_energy_discrepancies"=>delta[0]+2delta[2]+delta[4],
        "interpretation"=>"observed_method_agreement_not_rigorous_error_bounds",
        "chi_and_seed_dependence"=>"not_scanned")
    summary["magnetization_changes"] = [Dict(
        "from_Q"=>a,"to_Q"=>b,"delta_total_sz"=>(b-a)/2,
        "site_delta_sz"=>finals[b]["sz_profile"]-finals[a]["sz_profile"],
        "column_delta_sz"=>finals[b]["column_sz"]-finals[a]["column_sz"],
        "ed_column_delta_sz"=>finals[b]["projected_ed_column_sz"]-finals[a]["projected_ed_column_sz"],
        "edge_vs_bulk"=>"not_resolved_no_interior_column") for (a,b) in ((0,2),(2,4))]
    consistency = Dict{String,Any}("tolerance"=>1e-10,
        "ed_su2_energy_order"=>ed_energies[0] <= ed_energies[2]+1e-10 &&
                                  ed_energies[2] <= ed_energies[4]+1e-10,
        "dmrg_su2_energy_order"=>dmrg_energies[0] <= dmrg_energies[2]+1e-10 &&
                                    dmrg_energies[2] <= dmrg_energies[4]+1e-10,
        "max_magnetization_change_error"=>maximum(abs(sum(r["column_delta_sz"])-r["delta_total_sz"])
            for r in summary["magnetization_changes"]))
    consistency["passed"] = consistency["ed_su2_energy_order"] &&
        consistency["dmrg_su2_energy_order"] && consistency["max_magnetization_change_error"] <= 1e-10
    summary["consistency"] = consistency
end

function main(output)
    # The launcher owns execution.toml. Never reuse numerical output.
    isdir(output) || error("launch with run_static18.py and a new output directory")
    readdir(output) == ["execution.toml"] || error("output must contain only launcher execution.toml")
    execution = TOML.parsefile(joinpath(output,"execution.toml"))
    get(execution,"status","") == "running" || error("launcher is not running")
    compute_budget = execution["wall_limit_seconds"]-execution["termination_grace_seconds"]
    Threads.nthreads() == BLAS.get_num_threads() == 1 || error("one Julia/BLAS thread required")
    before = extra_hashes()
    lattice = kagome_cylinder(2,3)
    record = Dict{String,Any}("schema_version"=>1, "status"=>"running",
        "recorded_at_utc"=>string(now(UTC)), "scope"=>"static_N18_NN_sector_comparison_stage_A",
        "model"=>"nearest_neighbor_isotropic_Heisenberg_J1=1_J2=J3=hz=theta=0",
        "planned_Q"=>collect(CHARGES), "target_Q"=>2, "N"=>18,
        "code"=>KagomeDMRG._checkpoint_provenance(), "extra_source_sha256"=>before,
        "runtime"=>merge(KagomeDMRG._checkpoint_runtime(),Dict(
            "julia_threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads())),
        "solver"=>asdict(SOLVER), "ed_solver"=>asdict(ED_SOLVER),
        "ground_cluster_tolerance"=>CLUSTER_TOLERANCE, "limits"=>LIMITS,
        "retry_policy"=>"one_additional_six_sweeps_only_if_ED_comparison_fails",
        "time_limit_enforcement"=>"external_launcher_including_single_solves_and_startup",
        "wall_limit_seconds"=>execution["wall_limit_seconds"],
        "execution_record"=>"execution.toml", "variance_is_error_estimate"=>false,
        "bulk_plateau"=>"not_established", "quantized_pump"=>"not_measured",
        "phase_identification"=>"not_attempted", "runs"=>Dict{String,Any}[])
    persist() = KagomeDMRG._write_flux_record(joinpath(output,"validation.toml"),record)
    function phase!(run, phase)
        run["active_phase"] = phase
        record["worker_elapsed_seconds"] = elapsed()
        run["worker_remaining_seconds_at_phase_start"] = max(0.0,compute_budget-elapsed())
        persist()
        if phase in ("independent_ED","DMRG_6_sweeps","DMRG_additional_6_sweeps") &&
                elapsed() >= compute_budget
            throw(StudyTimeLimit())
        end
        println("Q=",run["Q"]," phase=",phase," elapsed=",round(elapsed();digits=2))
        flush(stdout)
    end
    persist()
    try
        for Q in CHARGES
            run = Dict{String,Any}("Q"=>Q,"status"=>"running",
                "configuration"=>KagomeDMRG._checkpoint_configuration(lattice,:seam,nothing;Q),
                "attempts"=>Dict{String,Any}[])
            push!(record["runs"],run)
            phase!(run,"independent_ED")
            ed_seconds = @elapsed ed = reference_eigensystem(2,3,0.0;Q,ED_SOLVER...)
            cluster = isempty(ed.values) ? Int[] :
                findall(e->e-first(ed.values)<CLUSTER_TOLERANCE,ed.values)
            boundary_resolved = !isempty(cluster) && length(cluster)<length(ed.values)
            reference = Dict{String,Any}("values"=>ed.values,
                "residual_norms"=>ed.residual_norms,"orthogonality_error"=>ed.orthogonality_error,
                "converged"=>ed.converged,"basis_dimension"=>length(ed.basis),
                "settings"=>asdict(ed.settings),"solver_info"=>asdict(ed.solver_info),
                "ground_cluster_size"=>length(cluster),"boundary_level_resolved"=>boundary_resolved,
                "ground_space_completeness"=>"partial_spectrum_not_proven",
                "elapsed_seconds"=>ed_seconds,"payload"=>save_reference(output,Q,ed))
            boundary_resolved && (reference["gap_above_ground_cluster"] =
                ed.values[length(cluster)+1]-first(ed.values))
            run["reference"] = reference
            run["reference_valid"] = ed.converged && boundary_resolved
            persist()
            if !run["reference_valid"]
                run["status"] = "unresolved_reference"
                phase!(run,"finished")
                continue
            end
            previous = nothing
            for attempt in 1:2
                phase!(run,attempt == 1 ? "DMRG_6_sweeps" : "DMRG_additional_6_sweeps")
                dmrg_seconds = @elapsed result = previous === nothing ?
                    run_dmrg(lattice,0.0;Q,SOLVER...) :
                    run_dmrg(lattice,0.0;Q,sites=previous.sites,psi0=previous.psi,SOLVER...)
                phase!(run,"diagnostics")
                diagnostics_seconds = @elapsed row = measure_state(result,ed)
                row["attempt"] = attempt
                row["total_sweeps"] = 6attempt
                row["dmrg_seconds"] = dmrg_seconds
                row["diagnostics_seconds"] = diagnostics_seconds
                push!(run["attempts"],row)
                phase!(run,"checkpoint")
                snapshot = save_checkpoint(output,result;baseline=result,theta_path=[0.0],status=:trial)
                loaded = load_checkpoint(snapshot,lattice;Q,status=:trial,expected_theta=0.0)
                row["checkpoint"] = relpath(snapshot,output)
                row["checkpoint_overlap_error"] = abs(abs(inner(loaded.psi,result.psi))-1)
                row["checkpoint_passed"] = loaded.Q == Q && row["checkpoint_overlap_error"] < 1e-12
                row["checkpoint_passed"] || error("checkpoint verification failed")
                persist()
                println("Q=",Q," sweeps=",row["total_sweeps"]," E=",row["energy"],
                    " residual=",row["residual_norm"]," passed=",row["passed"])
                flush(stdout)
                row["passed"] && break
                previous = result
            end
            run["status"] = last(run["attempts"])["passed"] ? "passed" : "unresolved_dmrg_accuracy"
            phase!(run,"finished")
            summarize!(record)
            persist()
            GC.gc()
        end
        summarize!(record)
        record["status"] = all(r["status"] == "passed" for r in record["runs"]) ?
            "passed_finite_cluster_comparison" : "completed_with_unresolved_sectors"
        if record["status"] == "passed_finite_cluster_comparison" &&
                !record["sector_comparison"]["consistency"]["passed"]
            record["status"] = "unresolved_sector_consistency"
        end
    catch exception
        record["status"] = exception isa StudyTimeLimit ? "time_budget_exceeded" :
            exception isa InterruptException ? "interrupted" : "failed"
        record["exception_type"] = string(nameof(typeof(exception)))
        for run in record["runs"]
            run["status"] == "running" && (run["status"] = record["status"])
        end
        rethrow()
    finally
        record["worker_elapsed_seconds"] = elapsed()
        record["extra_sources_unchanged"] = extra_hashes() == before
        record["extra_sources_unchanged"] || (record["status"] = "source_changed")
        persist()
    end
    println("Static N18 study: ",record["status"])
    return record["status"] == "passed_finite_cluster_comparison" ? 0 : 2
end

if abspath(PROGRAM_FILE) == @__FILE__
    length(ARGS) == 1 || error("usage: run_static18.py new-output-directory")
    exit(main(abspath(only(ARGS))))
end
