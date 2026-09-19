#!/usr/bin/env julia
# Reproduce same-runtime checkpoint restart in a separate Julia process.
# Run: julia --project=research --startup-file=no --threads=1 examples/validate_restart.jl
using SHA
using TOML
using Dates

const ROOT = normpath(joinpath(@__DIR__, ".."))
function active_environment_snapshot()
    project = Base.active_project()
    project === nothing && error("an explicit active project is required")
    project = abspath(project)
    manifest = Base.project_file_manifest_path(project)
    manifest !== nothing && isfile(project) && isfile(manifest) ||
        error("the active project needs an instantiated dependency manifest")
    manifest = abspath(manifest)
    return Dict{String,Any}("manifest"=>basename(manifest),
        "environment_sha256"=>Dict(
            "project:" * basename(project)=>bytes2hex(sha256(read(project))),
            "manifest:" * basename(manifest)=>bytes2hex(sha256(read(manifest)))))
end

const SOURCES = ["Project.toml", "research/Project.toml", "research/Manifest.toml",
    "research/Manifest-v1.13.toml",
    "src/KagomeDMRG.jl", "src/lattice.jl", "src/model.jl", "src/dmrg.jl",
    "src/observables.jl", "src/checkpoint.jl", "src/schmidt.jl", "src/continuation.jl",
    "test/reference_ed.jl", "test/itensor_helpers.jl", "examples/validate_restart.jl"]

function source_snapshot()
    hashes = Dict{String,String}()
    errors = String[]
    for path in SOURCES
        try
            hashes[path] = bytes2hex(sha256(read(joinpath(ROOT, path))))
        catch
            push!(errors, "source_hash_unavailable:" * path)
        end
    end
    revision = try
        strip(read(Cmd(Cmd(["git", "rev-parse", "HEAD"]); dir=ROOT), String))
    catch
        push!(errors, "git_revision_unavailable")
        "unavailable"
    end
    dirty = try
        !isempty(read(Cmd(Cmd(["git", "status", "--porcelain"]); dir=ROOT), String))
    catch
        push!(errors, "git_status_unavailable")
        "unavailable"
    end
    return merge(Dict("git_revision"=>revision, "worktree_dirty"=>dirty,
                "source_sha256"=>hashes, "errors"=>errors), active_environment_snapshot())
end
const SOURCE_BEFORE = source_snapshot()

using KagomeDMRG
using ITensors
using ITensorMPS
using KrylovKit
using LinearAlgebra

BLAS.set_num_threads(1)
ITensors.disable_threaded_blocksparse()

include(joinpath(ROOT, "test", "reference_ed.jl"))
include(joinpath(ROOT, "test", "itensor_helpers.jl"))

function write_record(path, record)
    ispath(path) && error("output already exists; choose a new output directory")
    mkpath(dirname(path))
    temporary, io = mktemp(dirname(path))
    try
        TOML.print(io, record; sorted=true)
        close(io)
        mv(temporary, path; force=false)
    finally
        isopen(io) && close(io)
        ispath(temporary) && rm(temporary)
    end
    return path
end

function point_metrics(result; bond=4)
    reference = reference_hamiltonian(1, 3, result.theta)
    eig = eigen(Hermitian(reference.H))
    selected = findall(e -> abs(e-eig.values[1]) < 1e-10, eig.values)
    ground = eig.vectors[:, selected]
    v = sector_amplitudes(result.psi, result.sites, reference.basis)
    projected = ground * (ground' * v)
    matched = projected / norm(projected)
    ed_sz = reference_sz(matched, reference.basis, 9)
    schmidt = schmidt_diagnostics(result.psi, bond)
    return Dict{String,Any}(
        "theta"=>result.theta, "energy"=>result.energy,
        "ed_energy"=>eig.values[1],
        "energy_error_per_site"=>abs(result.energy-eig.values[1])/9,
        "residual_norm"=>norm(reference.H*v-result.energy*v),
        "ground_space_leakage"=>norm(v-projected),
        "max_sz_error"=>maximum(abs.(result.sz-ed_sz)),
        "state_norm_error"=>abs(norm(v)-1),
        "total_sz_error"=>abs(sum(result.sz)-0.5),
        "absolute_variance"=>abs(result.variance),
        "schmidt_left_sz_error"=>abs(schmidt.mean_left_sz-sum(result.sz[1:bond])),
        "schmidt_probability_error"=>abs(sum(schmidt.probabilities)-1),
        "sz_profile"=>result.sz,
        "sweep_energies"=>result.sweep_energies,
        "measured_truncation_errors"=>result.max_truncation_errors,
        "solver"=>Dict(string(k)=>v for (k,v) in pairs(result.settings)),
        "schmidt"=>Dict("bond"=>bond, "region"=>"MPS_prefix_not_geometric_axis_cut",
            "probabilities"=>schmidt.probabilities, "left_q"=>schmidt.left_q,
            "entropy"=>schmidt.entropy, "mean_left_sz"=>schmidt.mean_left_sz))
end

function worker(checkpoint, output_directory, theta)
    lattice = kagome_cylinder(1, 3)
    loaded = load_checkpoint(checkpoint, lattice)
    resumed = resume_dmrg(checkpoint, lattice, theta;
                          expected_settings=loaded.settings)
    child_checkpoint = save_checkpoint(output_directory, resumed;
        baseline=loaded.baseline, theta_path=[loaded.theta_path; theta], status=:accepted)
    write_record(joinpath(output_directory, "worker.toml"), Dict(
        "checkpoint"=>relpath(child_checkpoint, output_directory),
        "theta"=>theta, "julia_threads"=>Threads.nthreads(),
        "blas_threads"=>BLAS.get_num_threads()))
end

function validate(output_directory)
    lattice = kagome_cylinder(1, 3)
    target = 2pi + 0.71
    limits = Dict("energy_error_per_site"=>1e-8, "residual_norm"=>2e-6,
        "ground_space_leakage"=>2e-6, "max_sz_error"=>1e-6,
        "state_norm_error"=>2e-12, "total_sz_error"=>2e-12,
        "absolute_variance"=>1e-10, "schmidt_left_sz_error"=>2e-12,
        "schmidt_probability_error"=>2e-12)
    restart_limits = Dict("energy_difference"=>1e-10, "max_sz_difference"=>1e-7,
        "overlap_magnitude_error"=>1e-10, "baseline_overlap_magnitude_error"=>2e-12,
        "schmidt_mean_left_sz_difference"=>1e-7)
    record = Dict{String,Any}(
        "schema_version"=>1, "recorded_at_utc"=>string(now(UTC)), "status"=>"error",
        "scope"=>"finite_size_checkpoint_and_schmidt_validation",
        "trajectory_status"=>"manual_fixture_not_adiabatic_branch_certification",
        "phase_identification"=>"not_attempted", "pump_quantization"=>"not_tested",
        "axial_pump"=>"not_testable_without_internal_geometric_cut_on_Lx_1",
        "model"=>Dict("name"=>"nearest_neighbor_kagome_heisenberg",
            "Jxy"=>1.0, "Jz"=>1.0, "longitudinal_field"=>0.0, "chirality_seed"=>0.0),
        "geometry"=>Dict("Lx"=>1, "Ly"=>3, "N"=>9, "bond_count"=>12,
            "wrap"=>[1.5, 3sqrt(3)/2], "ordering"=>"x_then_y_then_A_B_C",
            "termination"=>"complete_cells_with_outgoing_open_axis_bonds_removed",
            "axis_boundary"=>"open", "circumference_boundary"=>"periodic"),
        "conventions"=>Dict("Q"=>1, "physical_M"=>0.5, "Nup"=>5, "Ndown"=>4,
            "charge"=>"q=2Sz", "gauge"=>"seam",
            "exchange_phase"=>"S+_i_S-_j_has_exp(+im*wy*theta)"),
        "theta_path"=>[0.0, 0.37, target],
        "runtime"=>Dict("julia"=>string(VERSION), "julia_threads"=>Threads.nthreads(),
            "blas_threads"=>BLAS.get_num_threads(),
            "ITensors"=>string(Base.pkgversion(ITensors)),
            "ITensorMPS"=>string(Base.pkgversion(ITensorMPS)),
            "NDTensors"=>string(Base.pkgversion(ITensors.NDTensors)),
            "KrylovKit"=>string(Base.pkgversion(KrylovKit))),
        "acceptance_limits"=>limits, "restart_limits"=>restart_limits,
        "points"=>Dict{String,Any}[])
    function record_point!(label, result)
        metrics = point_metrics(result)
        metrics["label"] = label
        metrics["status"] = all(metrics[k] < v for (k,v) in limits) ? "passed" : "failed"
        push!(record["points"], metrics)
        println(label, " theta=", result.theta, " status=", metrics["status"],
                " energy_error/site=", metrics["energy_error_per_site"])
        return nothing
    end
    try
        baseline = run_dmrg(lattice, 0.0; seed=11)
        record_point!("baseline", baseline)
        point = run_dmrg(lattice, 0.37; sites=baseline.sites, psi0=baseline.psi, seed=11)
        record_point!("saved_point", point)
        accepted = save_checkpoint(joinpath(output_directory, "parent"), point;
            baseline, theta_path=[0.0, 0.37], status=:accepted)
        before = Dict(name=>bytes2hex(sha256(read(joinpath(accepted, name))))
                      for name in readdir(accepted))
        solver = (; (k=>v for (k,v) in pairs(point.settings) if k != :initialization)...)
        direct = run_dmrg(lattice, target; sites=point.sites, psi0=point.psi, solver...)
        record_point!("direct", direct)
        child_directory = joinpath(output_directory, "child")
        # All theta points are sequential. The child uses the same executable,
        # project, saved indices, solver settings, and one CPU/BLAS thread.
        command = `$(Base.julia_cmd()) --project=$(dirname(Base.active_project())) --startup-file=no --threads=1 $(@__FILE__) --worker $accepted $child_directory $target`
        run(command)
        worker_record = TOML.parsefile(joinpath(child_directory, "worker.toml"))
        restarted = load_checkpoint(joinpath(child_directory, worker_record["checkpoint"]),
                                    lattice; expected_theta=target,
                                    expected_settings=direct.settings)
        # Rebuild H independently of serialized state; no saved environments
        # or MPO objects are required to resume the solver.
        H = twisted_exchange_mpo(restarted.sites, lattice, target)
        energy = real(inner(restarted.psi', H, restarted.psi))
        variance = real(inner(H, restarted.psi, H, restarted.psi)) - energy^2
        resumed = (; psi=restarted.psi, sites=restarted.sites, theta=restarted.theta,
            energy, variance, sz=sz_profile(restarted.psi), settings=restarted.settings,
            sweep_energies=restarted.diagnostics.sweep_energies,
            max_truncation_errors=restarted.diagnostics.max_truncation_errors)
        record_point!("separate_process_restart", resumed)
        restart = Dict{String,Any}(
            "energy_difference"=>abs(direct.energy-resumed.energy),
            "max_sz_difference"=>maximum(abs.(direct.sz-resumed.sz)),
            "overlap_magnitude_error"=>abs(1-abs(inner(direct.psi, resumed.psi))),
            "baseline_overlap_magnitude_error"=>abs(1-abs(inner(restarted.baseline.psi,
                                                                baseline.psi))),
            "schmidt_mean_left_sz_difference"=>abs(
                schmidt_diagnostics(direct.psi, 4).mean_left_sz-
                schmidt_diagnostics(resumed.psi, 4).mean_left_sz),
            "sites_preserved"=>restarted.sites == point.sites,
            "charge_preserved"=>flux(restarted.psi) == QN("Sz", 1),
            "baseline_preserved"=>restarted.baseline.sz == baseline.sz,
            "theta_path_preserved"=>restarted.theta_path == [0.0, 0.37, target],
            "accepted_snapshot_unchanged"=>all(
                before[name] == bytes2hex(sha256(read(joinpath(accepted, name))))
                for name in keys(before)),
            "separate_process"=>true,
            "worker_julia_threads"=>worker_record["julia_threads"],
            "worker_blas_threads"=>worker_record["blas_threads"])
        exact_checks = ("sites_preserved", "charge_preserved", "baseline_preserved",
                        "theta_path_preserved", "accepted_snapshot_unchanged")
        restart["status"] = all(restart[k] < v for (k,v) in restart_limits) &&
            all(restart[k] for k in exact_checks) &&
            restart["worker_julia_threads"] == restart["worker_blas_threads"] == 1 ?
            "passed" : "failed"
        record["restart"] = restart
        record["status"] = all(p["status"] == "passed" for p in record["points"]) &&
            restart["status"] == "passed" ? "passed" : "failed"
    catch exception
        # Retain failure type, not arbitrary host/environment or process logs.
        record["exception_type"] = string(typeof(exception))
        showerror(stderr, exception)
        println(stderr)
    end
    source_after = source_snapshot()
    source_status = if !isempty(SOURCE_BEFORE["errors"]) || !isempty(source_after["errors"])
        "unavailable"
    elseif SOURCE_BEFORE["source_sha256"] != source_after["source_sha256"] ||
           SOURCE_BEFORE["environment_sha256"] != source_after["environment_sha256"] ||
           SOURCE_BEFORE["manifest"] != source_after["manifest"] ||
           SOURCE_BEFORE["git_revision"] != source_after["git_revision"]
        "source_changed"
    else
        "verified"
    end
    record["code"] = Dict("status"=>source_status,
                           "before"=>SOURCE_BEFORE, "after"=>source_after)
    source_status == "verified" || (record["status"] = "failed")
    output = write_record(joinpath(output_directory, "validation.toml"), record)
    println("Saved restart validation: ", output)
    record["status"] == "passed" || error("restart validation failed; outcomes were saved")
end

if !isempty(ARGS) && first(ARGS) == "--worker"
    length(ARGS) == 4 || error("worker requires checkpoint, output directory, and theta")
    worker(ARGS[2], ARGS[3], parse(Float64, ARGS[4]))
else
    length(ARGS) <= 1 || error("usage: validate_restart.jl [new-output-directory]")
    output_directory = if isempty(ARGS)
        parent = joinpath(ROOT, "outputs")
        mkpath(parent)
        mktempdir(parent; prefix="p1-restart-", cleanup=false)
    else
        path = abspath(only(ARGS))
        ispath(path) && error("output already exists; choose a new output directory")
        mkpath(path)
        path
    end
    validate(output_directory)
end
