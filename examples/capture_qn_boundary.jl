#!/usr/bin/env julia
# Isolated diagnostic process for the unmodified NDTensors 0.4.31 backend.
# Arguments: frozen KagomeDMRG source directory, new output directory.
# Use the frozen checkout's original dependency environment (historically its
# root --project), not the current research environment. No dependency is edited.
using ITensors, ITensorMPS, LinearAlgebra, Serialization, SHA, TOML, Dates
BLAS.set_num_threads(1)
ITensors.disable_threaded_blocksparse()

length(ARGS) == 2 || error("usage: capture_qn_boundary.jl frozen-source new-output-directory")
const SOURCE_ROOT = abspath(ARGS[1])
const OUTPUT_ROOT = abspath(ARGS[2])
ispath(OUTPUT_ROOT) && error("choose a new output directory")
mkpath(OUTPUT_ROOT)
include(joinpath(SOURCE_ROOT, "src", "KagomeDMRG.jl"))
using .KagomeDMRG

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

function source_hashes(root)
    files = ["Project.toml"]
    # Preserve either historical root or dedicated research lockfile layouts.
    append!(files, [path for path in ("Manifest.toml", "Manifest-v1.12.toml",
        "Manifest-v1.13.toml", "research/Project.toml", "research/Manifest.toml",
        "research/Manifest-v1.12.toml", "research/Manifest-v1.13.toml")
        if isfile(joinpath(root, path))])
    append!(files, [relpath(joinpath(dir, file), root)
        for (dir, _, names) in walkdir(joinpath(root, "src")) for file in names])
    return Dict(name => bytes2hex(sha256(read(joinpath(root, name)))) for name in files)
end

const BEFORE = source_hashes(SOURCE_ROOT)
const ENVIRONMENT_BEFORE = active_environment_snapshot()
const BACKEND_FILES = ("src/blocksparse/linearalgebra.jl",
    "src/lib/RankFactorization/src/truncate_spectrum.jl")
backend_hashes() = Dict(name => bytes2hex(sha256(read(joinpath(
    pkgdir(ITensors.NDTensors), name)))) for name in BACKEND_FILES)
const BACKEND_BEFORE = backend_hashes()
const CAPTURED = Ref(false)
const UPDATE = Ref(0)
const REPORT = Dict{String,Any}("status"=>"running", "schema_version"=>1,
    "scope"=>"unmodified_backend_failure_input_capture", "recorded_at_utc"=>string(now(UTC)),
    "source_sha256"=>BEFORE, "backend_sha256"=>BACKEND_BEFORE,
    "environment_before"=>ENVIRONMENT_BEFORE,
    "probe_sha256"=>bytes2hex(sha256(read(@__FILE__))),
    "runtime"=>Dict("julia"=>string(VERSION), "ITensors"=>string(pkgversion(ITensors)),
        "ITensorMPS"=>string(pkgversion(ITensorMPS)),
        "NDTensors"=>string(pkgversion(ITensors.NDTensors)),
        "julia_threads"=>Threads.nthreads(), "blas_threads"=>BLAS.get_num_threads()),
    "solver"=>Dict("seed"=>11, "nsweeps"=>6, "maxdim"=>256,
        "cutoff"=>0.0, "noise"=>0.0, "theta"=>0.0,
        "eigsolve_tol"=>1e-11, "eigsolve_krylovdim"=>40, "eigsolve_maxiter"=>20))

# A more-specific tracing method in this disposable process delegates unchanged
# to the upstream MPS update. The production package has no method replacement.
function ITensorMPS.replacebond!(PH::ITensorMPS.ProjMPO, psi::MPS,
                                bond::Int, phi::ITensor; kwargs...)
    UPDATE[] += 1
    before = deepcopy(phi)
    left = collect(commoninds(phi, psi[bond]))
    spec = ITensorMPS.replacebond!(psi, bond, phi; kwargs...)
    actual_rank = dim(linkind(psi, bond))
    reported_rank = length(eigs(spec))
    if actual_rank != reported_rank && !CAPTURED[]
        CAPTURED[] = true
        right = collect(uniqueinds(before, left))
        matrix = reshape(Array(before, left..., right...), prod(dim, left), :)
        dense_weights = svdvals(matrix).^2
        retained = Int(kwargs[:maxdim])
        expected_loss = sum(dense_weights[(retained+1):end]) / sum(dense_weights)
        reconstructed = psi[bond] * psi[bond+1]
        kept = reshape(Array(reconstructed, left..., right...), size(matrix))
        overlap_weight = abs2(dot(vec(matrix), vec(kept))) /
            (sum(abs2, matrix) * sum(abs2, kept))
        merge!(REPORT, Dict("bond"=>bond, "update"=>UPDATE[],
            "ortho"=>kwargs[:ortho], "actual_rank"=>actual_rank,
            "reported_rank"=>reported_rank, "reported_truncation_error"=>truncerror(spec),
            "input_norm_squared"=>sum(abs2, matrix),
            "dense_best_rank_loss"=>expected_loss,
            "normalized_projection_loss"=>max(0.0, 1-overlap_weight),
            "dense_weights"=>dense_weights,
            "spectrum"=>collect(eigs(spec)), "input_charge"=>val(flux(before), "Sz")))
        # Only the decomposition input and its explicit kwargs are saved.
        # No Hamiltonian environment, process state, or connection data is saved.
        payload = (; phi=before, left, kwargs=(; kwargs...), bond,
            spectrum=collect(eigs(spec)), actual_rank,
            reported_truncation_error=truncerror(spec))
        path = joinpath(OUTPUT_ROOT, "factorization_input.jls")
        open(io -> serialize(io, payload), path, "w")
        REPORT["payload_sha256"] = bytes2hex(sha256(read(path)))
        println("Captured truncation mismatch: bond=", bond,
            " reported_rank=", reported_rank, " retained_rank=", actual_rank)
    end
    return spec
end

try
    lattice = kagome_cylinder(2, 3)
    elapsed = @elapsed run_dmrg(lattice, 0.0; seed=11, nsweeps=6, maxdim=256,
        cutoff=0.0, noise=0.0, eigsolve_tol=1e-11, eigsolve_krylovdim=40,
        eigsolve_maxiter=20, outputlevel=1)
    REPORT["status"] = CAPTURED[] ? "unexpected_guard_bypass" : "failure_not_reproduced"
    REPORT["elapsed_seconds"] = elapsed
catch exception
    exception isa InterruptException && rethrow()
    exception isa OutOfMemoryError && rethrow()
    REPORT["status"] = CAPTURED[] ? "captured_guard_failure" : "uncaptured_failure"
    REPORT["exception_type"] = string(nameof(typeof(exception)))
end
REPORT["environment_after"] = active_environment_snapshot()
REPORT["sources_unchanged"] = source_hashes(SOURCE_ROOT) == BEFORE &&
    backend_hashes() == BACKEND_BEFORE && REPORT["environment_after"] == ENVIRONMENT_BEFORE
open(joinpath(OUTPUT_ROOT, "capture.toml"), "w") do io
    TOML.print(io, REPORT; sorted=true)
end
println("capture_status=", REPORT["status"])
REPORT["status"] == "captured_guard_failure" && REPORT["sources_unchanged"] ||
    error("expected boundary failure was not captured; diagnostic record retained")
