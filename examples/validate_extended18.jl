#!/usr/bin/env julia
# Two sequential small-system optimizer/observable checks, not a CSL pump.
# Run through run_extended18.py for an inclusive timeout inside every solve.
const WORKER_START = time_ns()
using KagomeDMRG, ITensors, ITensorMPS, LinearAlgebra, SHA, TOML, Dates, Serialization
BLAS.set_num_threads(1)
ITensors.disable_threaded_blocksparse()
const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "test", "reference_ed.jl"))
include(joinpath(ROOT, "test", "reference_extended_ed.jl"))
include(joinpath(ROOT, "test", "reference_eigensolve.jl"))
include(joinpath(ROOT, "test", "itensor_helpers.jl"))
include(joinpath(@__DIR__, "extended18_chirality_reference.jl"))

const COUPLINGS = (; J1=1.0, J2=0.5, J3=0.5)
const SOLVER = (; seed=11, nsweeps=6, maxdim=512, cutoff=0.0, noise=0.0,
    eigsolve_tol=1e-11, eigsolve_krylovdim=40, eigsolve_maxiter=20,
    measure_variance=true)
const ED_SOLVER = (; nev=4, blocksize=4, seed=5678, krylovdim=80,
    maxiter=200, tol=1e-11)
const CLUSTER_TOLERANCE = 1e-8
const LIMITS = Dict("energy_error_per_site"=>1e-8, "residual_norm"=>2e-6,
    "ground_space_leakage"=>2e-6, "max_sz_error"=>1e-6,
    "max_bond_energy_error"=>1e-6, "max_chirality_error"=>1e-6,
    "chirality_same_state_error"=>1e-10, "chirality_gauge_error"=>1e-10,
    "total_sz_error"=>1e-12, "state_norm_error"=>1e-12,
    "mpo_energy_error"=>1e-10, "bond_energy_sum_error"=>1e-10,
    "schmidt_density_error"=>1e-10, "transfer_conservation_error"=>1e-10,
    "schmidt_transfer_error"=>1e-10, "ed_transfer_error"=>1e-6,
    "variance_residual_squared_error"=>1e-9)
const EXTRA_SOURCES = ("examples/validate_extended18.jl", "examples/run_extended18.py",
    "examples/run_static18.py", "examples/run_static27.py",
    "examples/extended18_chirality_reference.jl", "test/reference_ed.jl",
    "test/reference_extended_ed.jl", "test/reference_eigensolve.jl", "test/itensor_helpers.jl")
extra_hashes() = Dict(p=>bytes2hex(sha256(read(joinpath(ROOT,p)))) for p in EXTRA_SOURCES)
asdict(value) = Dict(string(k)=>v for (k,v) in pairs(value))
elapsed() = (time_ns()-WORKER_START)/1e9
bond_key(b) = b.i < b.j ? (b.i,b.j,b.wy) : (b.j,b.i,-b.wy)

function save_reference(output, point, theta, ed)
    path = joinpath(output,"reference_point$(point).jls")
    temporary, io = mktemp(output)
    try
        serialize(io,(; Q=0,theta,basis=ed.basis,values=ed.values,vectors=ed.vectors))
        close(io)
        mv(temporary,path)
    finally
        isopen(io) && close(io)
        isfile(temporary) && rm(temporary)
    end
    return Dict("file"=>basename(path),"sha256"=>bytes2hex(sha256(read(path))))
end

function measure_state(result, ed, baseline, baseline_row)
    lattice, theta = result.lattice, result.theta
    v = sector_amplitudes(result.psi,result.sites,ed.basis)
    cluster = findall(e->e-first(ed.values)<CLUSTER_TOLERANCE,ed.values)
    ground = ed.vectors[:,cluster]
    projection = ground*(ground'*v)
    norm(projection)>0 || error("zero projection onto reference ground cluster")
    matched = projection/norm(projection)
    ed_profile = reference_sz(matched,ed.basis,18)
    action = ed.H*v
    rayleigh = real(dot(v,action))/real(dot(v,v))
    residual = norm(action-result.energy*v)
    bond_seconds = @elapsed energies = bond_energies(result.psi,lattice,theta)
    independent_bonds = reference_extended_bonds(2,3;COUPLINGS...)
    expected_bonds = Dict(bond_key(b)=>
        b.Jz*real(reference_correlation(matched,ed.basis,18,b.i,b.j)) +
        b.Jxy*real(cis(b.wy*theta)*reference_correlation(matched,ed.basis,18,b.i,b.j;
            operators=("S+","S-"))) for b in independent_bonds)
    schmidt_seconds = @elapsed schmidt = schmidt_diagnostics(result.psi,9)
    base = baseline === nothing ? result : baseline
    base_ed_sz = baseline_row === nothing ? ed_profile : baseline_row["projected_ed_sz_profile"]
    base_mean = baseline_row === nothing ? schmidt.mean_left_sz : baseline_row["schmidt_mean_left_sz"]
    transfer = only(spin_transfer(lattice,base.sz,result.sz))
    ed_transfer = sum((ed_profile-base_ed_sz)[10:18])
    triangles = oriented_triangles(lattice)
    chirality_seconds = @elapsed chiralities = triangle_chiralities(result.psi,lattice,theta)
    reference_chirality_seconds = @elapsed begin
        same_chirality = extended18_reference_chiralities(v,ed.basis,18,triangles,theta)
        matched_chirality = extended18_reference_chiralities(matched,ed.basis,18,triangles,theta)
    end
    uniform_psi = deepcopy(result.psi)
    for (i,chi) in enumerate(gauge_angles(lattice,theta))
        s = result.sites[i]
        U = cos(chi/2)*op("Id",s) + 2im*sin(chi/2)*op("Sz",s)
        uniform_psi[i] = noprime(U*uniform_psi[i])
    end
    gauge_seconds = @elapsed uniform_chirality =
        triangle_chiralities(uniform_psi,lattice,theta;gauge=:uniform)
    families = bond_families(lattice)
    row = Dict{String,Any}(
        "energy"=>result.energy,"ed_energy"=>first(ed.values),
        "energy_error_per_site"=>abs(result.energy-first(ed.values))/18,
        "residual_norm"=>residual,"ground_space_leakage"=>norm(v-projection),
        "ground_projection_norm"=>norm(projection),
        "max_sz_error"=>maximum(abs.(result.sz-ed_profile)),
        "max_bond_energy_error"=>maximum(abs(energies[k]-expected_bonds[bond_key(b)])
            for (k,b) in enumerate(lattice.bonds)),
        "max_chirality_error"=>maximum(abs.(chiralities-matched_chirality)),
        "chirality_same_state_error"=>maximum(abs.(chiralities-same_chirality)),
        "chirality_gauge_error"=>maximum(abs.(chiralities-uniform_chirality)),
        "total_sz_error"=>abs(sum(result.sz)),"state_norm_error"=>abs(norm(v)-1),
        "independent_rayleigh_energy"=>rayleigh,
        "independent_rayleigh_residual_norm"=>norm(action-rayleigh*v),
        "mpo_energy_error"=>abs(rayleigh-result.energy),
        "sz_profile"=>result.sz,"projected_ed_sz_profile"=>ed_profile,
        "column_sz"=>[sum(result.sz[1:9]),sum(result.sz[10:18])],
        "bond_energy_profile"=>energies,
        "family_energy"=>Dict(string(f)=>sum(energies[families .== f]) for f in (:J1,:J2,:J3)),
        "bond_energy_sum_error"=>abs(sum(energies)-result.energy),
        "chirality_profile"=>chiralities,"projected_ed_chirality_profile"=>matched_chirality,
        "same_state_reference_chirality_profile"=>same_chirality,
        "chirality_mean"=>sum(chiralities)/length(chiralities),
        "raw_variance"=>result.variance,
        "variance_residual_squared_error"=>abs(result.variance-residual^2),
        "variance_roundoff_scale"=>100eps(Float64)*max(1,result.energy^2),
        "maxlinkdim"=>maxlinkdim(result.psi),"sweep_energies"=>result.sweep_energies,
        "last_sweep_energy_change"=>abs(result.sweep_energies[end]-result.sweep_energies[end-1]),
        "measured_truncation_errors"=>result.max_truncation_errors,
        "schmidt_entropy"=>schmidt.entropy,"schmidt_mean_left_sz"=>schmidt.mean_left_sz,
        "schmidt_variance_left_sz"=>schmidt.variance_left_sz,
        "schmidt_charge_weights"=>[Dict("left_q"=>q,
            "probability"=>sum(schmidt.probabilities[schmidt.left_q .== q]))
            for q in sort(unique(schmidt.left_q))],
        "schmidt_density_error"=>abs(schmidt.mean_left_sz-sum(result.sz[1:9])),
        "transfer_left"=>transfer.left,"transfer_right"=>transfer.right,
        "transfer_conservation_error"=>abs(transfer.total),
        "schmidt_transfer_right"=>-(schmidt.mean_left_sz-base_mean),
        "schmidt_transfer_error"=>abs(transfer.right+schmidt.mean_left_sz-base_mean),
        "projected_ed_transfer_right"=>ed_transfer,"ed_transfer_error"=>abs(transfer.right-ed_transfer),
        "overlap_with_zero_flux"=>abs(inner(base.psi,result.psi))/(norm(base.psi)*norm(result.psi)),
        "measurement_seconds"=>Dict("bond_energies"=>bond_seconds,"schmidt"=>schmidt_seconds,
            "chirality"=>chirality_seconds,"reference_chirality"=>reference_chirality_seconds,
            "uniform_chirality"=>gauge_seconds),"settings"=>asdict(result.settings))
    row["failed_metrics"] = sort([k for (k,limit) in LIMITS if !isfinite(row[k]) || row[k]>limit])
    row["passed"] = isempty(row["failed_metrics"])
    return row
end

function main(output)
    isdir(output) && readdir(output)==["execution.toml"] ||
        error("launch with run_extended18.py and a new output directory")
    execution = TOML.parsefile(joinpath(output,"execution.toml"))
    get(execution,"status","")=="running" || error("launcher is not running")
    Threads.nthreads()==BLAS.get_num_threads()==1 || error("one Julia/BLAS thread required")
    before = extra_hashes()
    lattice = kagome_j1j2j3_cylinder(2,3;COUPLINGS...)
    record = Dict{String,Any}("schema_version"=>1,"status"=>"running",
        "recorded_at_utc"=>string(now(UTC)),"scope"=>"extended_N18_Q0_ED_DMRG_measurement_validation",
        "model"=>"J1_J2_J3_hexagon_opposite_only_no_field_or_chirality_seed",
        "couplings"=>asdict(COUPLINGS),"N"=>18,"Q"=>0,"planned_theta"=>[0.0,0.37],
        "configuration"=>KagomeDMRG._checkpoint_configuration(lattice,:seam,nothing;Q=0),
        "code"=>KagomeDMRG._checkpoint_provenance(),"extra_source_sha256"=>before,
        "runtime"=>merge(KagomeDMRG._checkpoint_runtime(),Dict(
            "julia_threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads())),
        "solver"=>asdict(SOLVER),"ed_solver"=>asdict(ED_SOLVER),"limits"=>LIMITS,
        "ground_cluster_tolerance"=>CLUSTER_TOLERANCE,
        "retry_policy"=>"one_additional_six_sweeps_at_same_theta_if_ED_comparison_fails",
        "time_limit_enforcement"=>"external_launcher_including_single_solves_and_startup",
        "wall_limit_seconds"=>execution["wall_limit_seconds"],"execution_record"=>"execution.toml",
        "quantized_pump"=>"not_measured","trajectory_acceptance"=>"not_attempted_trial_points_only",
        "bulk_region"=>"none_two_columns_one_axial_cut","phase_identification"=>"not_attempted",
        "triangles"=>[Dict("sites"=>collect(t.sites),"image_y"=>collect(t.image_y),
            "kind"=>string(t.kind)) for t in oriented_triangles(lattice)],
        "runs"=>Dict{String,Any}[])
    persist() = KagomeDMRG._write_flux_record(joinpath(output,"validation.toml"),record)
    function phase!(run,phase)
        run["active_phase"] = phase
        record["worker_elapsed_seconds"] = elapsed()
        persist()
        println("theta=",run["theta"]," phase=",phase," elapsed=",round(elapsed();digits=2))
        flush(stdout)
    end
    baseline = baseline_row = nothing
    persist()
    try
        selfcheck_seconds = @elapsed selfcheck = extended18_reference_chirality_selfcheck()
        record["chirality_reference_selfcheck"] = asdict(selfcheck)
        record["chirality_reference_selfcheck_seconds"] = selfcheck_seconds
        persist()
        for (point,theta) in enumerate(record["planned_theta"])
            run = Dict{String,Any}("theta"=>theta,"status"=>"running","attempts"=>Dict{String,Any}[])
            push!(record["runs"],run)
            phase!(run,"independent_ED")
            ed_seconds = @elapsed begin
                ref = reference_extended_hamiltonian(2,3,theta;Q=0,sparse_matrix=true,COUPLINGS...)
                low = _reference_low_eigensystem(ref.H;ED_SOLVER...)
                ed = (;ref...,low...)
            end
            cluster = findall(e->e-first(ed.values)<CLUSTER_TOLERANCE,ed.values)
            resolved = !isempty(cluster) && length(cluster)<length(ed.values)
            run["reference"] = Dict("values"=>ed.values,"residual_norms"=>ed.residual_norms,
                "orthogonality_error"=>ed.orthogonality_error,"converged"=>ed.converged,
                "basis_dimension"=>length(ed.basis),"settings"=>asdict(ed.settings),
                "solver_info"=>asdict(ed.solver_info),"observed_ground_cluster_size"=>length(cluster),
                "boundary_level_resolved"=>resolved,"ground_space_completeness"=>"partial_spectrum_not_proven",
                "elapsed_seconds"=>ed_seconds,"payload"=>save_reference(output,point,theta,ed))
            resolved && (run["reference"]["finite_cluster_level_spacing"] =
                ed.values[length(cluster)+1]-first(ed.values))
            if !(ed.converged && resolved)
                run["status"] = "unresolved_reference"
                record["status"] = "unresolved_reference"
                return record
            end
            previous = baseline
            for attempt in 1:2
                phase!(run,attempt==1 ? "DMRG_6_sweeps" : "DMRG_additional_6_sweeps")
                function progress(event)
                    if event.kind==:sweep
                        run["last_sweep"] = asdict(event)
                        run["last_sweep"]["kind"] = string(event.kind)
                        run["last_sweep"]["attempt"] = attempt
                        persist()
                        println("theta=",theta," attempt=",attempt," sweep=",event.sweep,
                            " energy=",event.energy," maxdim=",event.maxlinkdim)
                        flush(stdout)
                    end
                end
                dmrg_seconds = @elapsed result = previous===nothing ?
                    run_dmrg(lattice,theta;Q=0,progress_callback=progress,SOLVER...) :
                    run_dmrg(lattice,theta;Q=0,sites=previous.sites,psi0=previous.psi,
                        progress_callback=progress,SOLVER...)
                phase!(run,"checkpoint")
                snapshot = save_checkpoint(output,result;baseline=baseline===nothing ? result : baseline,
                    theta_path=point==1 ? [0.0] : [0.0,theta],status=:trial)
                run["latest_checkpoint"] = relpath(snapshot,output)
                persist()
                phase!(run,"diagnostics")
                diagnostics_seconds = @elapsed row = measure_state(result,ed,baseline,baseline_row)
                row["attempt"] = attempt
                row["total_sweeps_at_theta"] = 6attempt
                row["dmrg_seconds"] = dmrg_seconds
                row["diagnostics_seconds"] = diagnostics_seconds
                row["checkpoint"] = relpath(snapshot,output)
                push!(run["attempts"],row)
                phase!(run,"reload_verification")
                loaded = load_checkpoint(snapshot,lattice;Q=0,status=:trial,expected_theta=theta)
                row["checkpoint_overlap_error"] = abs(abs(inner(loaded.psi,result.psi))-1)
                row["checkpoint_baseline_profile_error"] = maximum(abs.(loaded.baseline.sz-
                    (baseline===nothing ? result.sz : baseline.sz)))
                row["checkpoint_passed"] = loaded.Q==0 && row["checkpoint_overlap_error"]<1e-12 &&
                    row["checkpoint_baseline_profile_error"]<1e-12
                row["checkpoint_passed"] || error("checkpoint reload verification failed")
                persist()
                println("theta=",theta," attempt=",attempt," residual=",row["residual_norm"],
                    " chirality_error=",row["chirality_same_state_error"]," passed=",row["passed"])
                flush(stdout)
                if row["passed"]
                    if point==1
                        baseline,baseline_row = result,row
                    end
                    break
                end
                previous = result
            end
            run["status"] = last(run["attempts"])["passed"] ? "passed" : "outside_limits"
            phase!(run,"finished")
            if run["status"]!="passed"
                record["status"] = "outside_limits"
                return record
            end
        end
        record["status"] = "passed"
    catch exception
        record["status"] = exception isa InterruptException ? "interrupted" : "failed"
        record["exception_type"] = string(nameof(typeof(exception)))
        rethrow()
    finally
        record["worker_elapsed_seconds"] = elapsed()
        record["extra_sources_unchanged"] = extra_hashes()==before
        record["extra_sources_unchanged"] || (record["status"]="source_changed")
        persist()
    end
    return record
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==1 || error("usage: validate_extended18.jl new-output-directory")
    record = main(abspath(only(ARGS)))
    record["status"]=="passed" || error("extended18 study unresolved; inspect validation.toml")
end
