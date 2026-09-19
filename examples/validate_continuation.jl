#!/usr/bin/env julia
# Bounded P2 driver validation, including an intentionally unresolved branch.
# julia --project=research --startup-file=no --threads=1 examples/validate_continuation.jl
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
    "src/observables.jl", "src/schmidt.jl", "src/checkpoint.jl", "src/continuation.jl",
    "test/reference_ed.jl", "test/itensor_helpers.jl", "examples/validate_continuation.jl"]
function source_hashes()
    hashes = Dict(path=>bytes2hex(sha256(read(joinpath(ROOT, path)))) for path in SOURCES)
    return merge(Dict("source_sha256"=>hashes), active_environment_snapshot())
end
const SOURCE_BEFORE = source_hashes()

using KagomeDMRG
using ITensors
using ITensorMPS
using LinearAlgebra

BLAS.set_num_threads(1)
ITensors.disable_threaded_blocksparse()
include(joinpath(ROOT, "test", "reference_ed.jl"))
include(joinpath(ROOT, "test", "itensor_helpers.jl"))

function policy(; kwargs...)
    limits = (; min_overlap=0.0, max_density_change=1.0, max_entropy_change=10.0,
        max_schmidt_change=9.0, max_variance=1e-9, max_truncation_error=1e-10,
        max_sweep_energy_change=1e-9, max_cut_spread=1e-9, consistency_tol=1e-9)
    return FluxPolicy(; merge(limits, (; kwargs...))...)
end

function ed_metrics(saved)
    reference = reference_hamiltonian(1, 3, saved.theta)
    eig = eigen(Hermitian(reference.H))
    ground = eig.vectors[:, findall(e -> abs(e-eig.values[1]) < 1e-9, eig.values)]
    v = sector_amplitudes(saved.psi, saved.sites, reference.basis)
    projected = ground * (ground' * v)
    matched = projected / norm(projected)
    return Dict("theta"=>saved.theta,
        "energy_error_per_site"=>abs(saved.diagnostics.energy-eig.values[1])/9,
        "residual_norm"=>norm(reference.H*v-saved.diagnostics.energy*v),
        "ground_space_leakage"=>norm(v-projected),
        "max_sz_error"=>maximum(abs.(saved.diagnostics.sz-reference_sz(matched,reference.basis,9))),
        "ground_degeneracy"=>size(ground,2))
end

function ed_points(trajectory, lattice; rejected=false)
    points = [ed_metrics(load_checkpoint(path, lattice)) for path in trajectory.accepted_checkpoints]
    if rejected
        for trial in trajectory.trials
            isempty(trial["checkpoint"]) && continue
            path = joinpath(dirname(trajectory.output_path), trial["checkpoint"])
            push!(points, ed_metrics(load_checkpoint(path, lattice; status=:trial)))
        end
    end
    return points
end

function trajectory_record(trajectory, directory)
    journal = TOML.parsefile(trajectory.output_path)
    return Dict{String,Any}("status"=>string(trajectory.status), "reason"=>trajectory.reason,
        "journal"=>relpath(trajectory.output_path,directory),
        "theta_path"=>trajectory.theta_path, "configuration"=>journal["configuration"],
        "policy"=>journal["policy"], "solver"=>journal["solver"],
        "initial_step"=>journal["initial_step"], "min_step"=>journal["min_step"],
        "cuts"=>journal["cuts"], "diagnostic_bonds"=>journal["diagnostic_bonds"],
        "bulk_sites"=>journal["bulk_sites"], "axial_cut_coverage"=>journal["axial_cut_coverage"],
        "initial_diagnostics"=>journal["initial"]["diagnostics"],
        "trials"=>[Dict("theta"=>t["theta"], "from_theta"=>t["from_theta"],
            "step"=>t["step"], "status"=>t["status"], "reasons"=>t["reasons"],
            "diagnostics"=>t["diagnostics"]) for t in trajectory.trials])
end

function real_zero_state(lattice)
    # This real ED vector is an arbitrary state in the degenerate zero-flux
    # ground space, used to expose a physical obstruction to overlap refinement.
    reference = reference_hamiltonian(1,3,0.0)
    eig = eigen(Hermitian(reference.H))
    amplitudes = real.(eig.vectors[:,1])
    amplitudes ./= norm(amplitudes)
    sites = spin_sites(lattice)
    tensor = ITensor(ComplexF64, sites...)
    for (amplitude, bits) in zip(amplitudes, reference.basis)
        tensor[[sites[j] => (isodd(bits >> (j-1)) ? 1 : 2) for j in 1:9]...] = amplitude
    end
    return run_dmrg(lattice,0.0; sites,psi0=MPS(tensor,sites),seed=11)
end

function validate(directory)
    record = Dict{String,Any}("schema_version"=>1,"status"=>"error",
        "recorded_at_utc"=>string(now(UTC)),
        "scope"=>"adaptive_continuation_controls_not_phase_identification",
        "pump_quantization"=>"not_tested", "positive_CSL_control"=>"not_implemented",
        "runtime"=>merge(KagomeDMRG._checkpoint_runtime(), Dict(
            "julia_threads"=>Threads.nthreads(), "blas_threads"=>BLAS.get_num_threads())),
        "code"=>Dict{String,Any}("before"=>SOURCE_BEFORE,
            "provenance"=>KagomeDMRG._checkpoint_provenance()),
        "runs"=>Dict{String,Any}())
    try
        # Full axial geometry with two cuts and an analytically unique product
        # state. Jxy=Jz=0, H=-sum(hz[i]*Sz[i]); this is an altered-model control.
        lattice = kagome_cylinder(3,3;Jxy=0,Jz=0)
        sites = spin_sites(lattice)
        labels = [fill("Up",15);fill("Dn",12)]
        hz = [fill(1.0,15);fill(-1.0,12)]
        baseline = run_dmrg(lattice,0.0;sites,psi0=MPS(ComplexF64,sites,labels),hz,
            nsweeps=2,maxdim=2,cutoff=0.0)
        zero_policy = policy(min_overlap=1-1e-10,max_density_change=1e-10,
            max_entropy_change=1e-10,max_schmidt_change=1e-10)
        first = continue_flux(lattice,[2pi];start=baseline,output_root=directory,
            policy=zero_policy,initial_step=2pi,min_step=pi/32,hz,
            cuts=[1,2],diagnostic_bonds=[9,18],bulk_sites=10:18)
        rest = continue_flux(lattice,[4pi,6pi,0.0];start=first.last_checkpoint,
            output_root=directory,policy=zero_policy,initial_step=2pi,min_step=pi/32,hz,
            cuts=[1,2],diagnostic_bonds=[9,18],bulk_sites=10:18)
        zero = trajectory_record(rest,directory)
        zero["first_leg"] = trajectory_record(first,directory)
        zero["max_abs_transfer"] = maximum(abs(t[key]) for run in (first,rest)
            for point in run.trials for t in point["diagnostics"]["transfers"]
            for key in ("left","right","schmidt_left","schmidt_right"))
        zero["max_cut_spread"] = maximum(p["diagnostics"]["cut_spread"] for run in (first,rest) for p in run.trials)
        zero["max_energy_error"] = maximum(abs(load_checkpoint(path,lattice;hz).diagnostics.energy+13.5)
            for run in (first,rest) for path in run.accepted_checkpoints)
        zero["expected_energy"] = -13.5
        zero["status_check"] = first.status == rest.status == :completed &&
            rest.theta_path == [0.0,2pi,4pi,6pi,4pi,2pi,0.0] &&
            zero["max_abs_transfer"] < 1e-12 && zero["max_energy_error"] < 1e-12
        record["runs"]["longitudinal_zero_control"] = zero
        println("zero control: ",rest.status," max transfer=",zero["max_abs_transfer"])

        small = kagome_cylinder(1,3)
        baseline9 = run_dmrg(small,0.0;seed=11)
        algebra = continue_flux(small,[0.17,0.37];start=baseline9,output_root=directory,
            policy=policy(),initial_step=0.2,min_step=0.0125,diagnostic_bonds=[4])
        algebra_record = trajectory_record(algebra,directory)
        algebra_record["interpretation"] = "Permissive overlap validates driver algebra and ED agreement only; no axial cut or adiabatic branch claim."
        algebra_record["ed_points"] = ed_points(algebra,small)
        algebra_record["status_check"] = algebra.status == :completed &&
            all(p["energy_error_per_site"] < 1e-8 && p["residual_norm"] < 2e-6 &&
                p["ground_space_leakage"] < 2e-6 && p["max_sz_error"] < 1e-6
                for p in algebra_record["ed_points"])
        record["runs"]["nine_site_ed"] = algebra_record
        println("nine-site ED: ",algebra.status)

        degenerate = real_zero_state(small)
        unresolved = continue_flux(small,[0.01];start=degenerate,output_root=directory,
            policy=policy(min_overlap=0.95),initial_step=0.01,min_step=0.0025,
            diagnostic_bonds=[4])
        negative = trajectory_record(unresolved,directory)
        negative["interpretation"] = "Real zero-flux ground state in a twofold degenerate space cannot meet the positive-flux overlap limit by step refinement alone. Valid rejected states are retained; this is not topological evidence."
        negative["ed_points"] = ed_points(unresolved,small;rejected=true)
        negative["status_check"] = unresolved.status == :unresolved &&
            unresolved.reason == "min_step" && unresolved.theta_path == [0.0] &&
            length(unresolved.trials) == 3 && all("low_overlap" in p["reasons"] for p in unresolved.trials) &&
            all(p["energy_error_per_site"] < 1e-8 && p["residual_norm"] < 2e-6
                for p in negative["ed_points"])
        record["runs"]["degenerate_initial_state"] = negative
        println("degenerate initial state: ",unresolved.status," reason=",unresolved.reason)
        record["status"] = all(r["status_check"] for r in values(record["runs"])) ? "passed" : "failed"
    catch exception
        record["exception_type"] = string(nameof(typeof(exception)))
        showerror(stderr,exception)
        println(stderr)
    end
    after = source_hashes()
    record["code"]["after"] = after
    record["code"]["status"] = after == SOURCE_BEFORE ? "verified" : "source_changed"
    after == SOURCE_BEFORE || (record["status"] = "failed")
    path = joinpath(directory,"validation.toml")
    ispath(path) && error("validation output already exists")
    KagomeDMRG._write_flux_record(path,record)
    println("Saved continuation validation: ",path)
    record["status"] == "passed" || error("validation failed; outcomes retained")
end

length(ARGS) <= 1 || error("usage: validate_continuation.jl [new-output-directory]")
directory = if isempty(ARGS)
    parent = joinpath(ROOT,"outputs")
    mkpath(parent)
    mktempdir(parent;prefix="p2-continuation-",cleanup=false)
else
    path = abspath(only(ARGS))
    ispath(path) && error("output already exists; choose a new directory")
    mkpath(path)
    path
end
validate(directory)
