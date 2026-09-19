#!/usr/bin/env julia
# Reproduce the small-system numerical checks with an allowlisted TOML record.
# Run from the repository: julia --project=research --startup-file=no examples/validate_small_system.jl
using SHA
using TOML
using Dates

const ROOT = normpath(joinpath(@__DIR__, ".."))
# Take the source snapshot before loading project code, and compare it again
# after numerics. Never attribute a running calculation to later edits.
# Track both research baselines; record the active environment separately.
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
    "src/KagomeDMRG.jl", "src/lattice.jl",
    "src/model.jl", "src/dmrg.jl", "src/observables.jl", "src/checkpoint.jl",
    "src/schmidt.jl", "src/continuation.jl", "test/reference_ed.jl",
    "test/itensor_helpers.jl", "examples/validate_small_system.jl"]
function provenance_snapshot()
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
const SOURCE_BEFORE = provenance_snapshot()

using KagomeDMRG
using ITensors
using ITensorMPS
using KrylovKit
using LinearAlgebra

BLAS.set_num_threads(1)
ITensors.disable_threaded_blocksparse()

include(joinpath(ROOT, "test", "reference_ed.jl"))
include(joinpath(ROOT, "test", "itensor_helpers.jl"))

length(ARGS) <= 1 || error("usage: validate_small_system.jl [new-output.toml]")
output = if isempty(ARGS)
    directory = joinpath(ROOT, "outputs")
    mkpath(directory)
    joinpath(mktempdir(directory; prefix="p0-reference-", cleanup=false), "validation.toml")
else
    abspath(only(ARGS))
end
ispath(output) && error("output already exists; choose a new path to preserve prior runs")
mkpath(dirname(output))

limits = Dict("energy_error_per_site" => 1e-8, "residual_norm" => 2e-6,
    "ground_space_leakage" => 2e-6, "max_sz_error" => 1e-6,
    "max_zz_error" => 1e-6, "max_pm_error" => 1e-6,
    "state_norm_error" => 1e-12, "absolute_variance" => 1e-10,
    "total_sz_error" => 1e-12)
points = Dict{String,Any}[]
for (theta, seed) in ((0.0, 11), (0.37, 11), (0.37, 29))
    try
        check = check_small_system(theta; seed)
        m, result = check.metrics, check.result
        point = Dict{String,Any}(string(k) => v for (k,v) in pairs(m))
        passed = all(point[name] < limits[name] for name in
            ("energy_error_per_site", "residual_norm", "ground_space_leakage",
             "max_sz_error", "max_zz_error", "max_pm_error", "state_norm_error")) &&
            abs(m.variance) < limits["absolute_variance"] &&
            abs(m.total_sz-0.5) < limits["total_sz_error"] &&
            m.ground_degeneracy == (theta == 0 ? 2 : 1)
        point["status"] = passed ? "passed" : "failed"
        point["sz_profile"] = result.sz
        point["sweep_energies"] = result.sweep_energies
        point["measured_truncation_errors"] = result.max_truncation_errors
        point["solver"] = Dict(string(k) => v for (k,v) in pairs(result.settings))
        push!(points, point)
        println("theta=", theta, " seed=", seed, " status=", point["status"],
                " energy_error/site=", m.energy_error_per_site,
                " max_Sz_error=", m.max_sz_error)
    catch exception
        # Keep failed points; do not archive arbitrary process/environment logs.
        push!(points, Dict("theta"=>theta, "seed"=>seed, "status"=>"error",
                           "exception_type"=>string(typeof(exception))))
        showerror(stderr, exception)
        println(stderr)
    end
end

# Record only explicit reproducibility fields, never environment variables or
# connection settings. A provenance error must not discard computed outcomes.
source_after = provenance_snapshot()
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
record = Dict{String,Any}(
    "schema_version"=>1, "recorded_at_utc"=>string(now(UTC)),
    "status"=>((source_status == "verified" && all(p["status"] == "passed" for p in points)) ? "passed" : "failed"),
    "scope"=>"finite_size_reference_backend_validation",
    "point_mode"=>"independent_single_points_not_flux_continuation",
    "ground_state_comparison"=>"project_DMRG_vector_into_full_ED_ground_space",
    "phase_identification"=>"not_attempted", "pump_quantization"=>"not_tested",
    "model"=>Dict("name"=>"nearest_neighbor_kagome_heisenberg", "Jxy"=>1.0,
                  "Jz"=>1.0, "longitudinal_field"=>0.0, "chirality_seed"=>0.0),
    "geometry"=>Dict("Lx"=>1, "Ly"=>3, "N"=>9, "bond_count"=>12,
        "a1"=>[1.0,0.0], "a2"=>[0.5,sqrt(3)/2], "wrap"=>[1.5,3sqrt(3)/2],
        "axis_boundary"=>"open", "circumference_boundary"=>"periodic",
        "termination"=>"complete_cells_with_outgoing_open_axis_bonds_removed",
        "ordering"=>"x_then_y_then_A_B_C", "gauge"=>"seam",
        "exchange_phase"=>"S+_i_S-_j_has_exp(+im*wy*theta)"),
    "sector"=>Dict("Q"=>1, "physical_M"=>0.5, "Nup"=>5, "Ndown"=>4,
                    "charge_convention"=>"q=2Sz"),
    "runtime"=>Dict("julia"=>string(VERSION), "julia_threads"=>Threads.nthreads(),
                     "blas_threads"=>BLAS.get_num_threads(),
                     "ITensors"=>string(Base.pkgversion(ITensors)),
                     "ITensorMPS"=>string(Base.pkgversion(ITensorMPS)),
                     "NDTensors"=>string(Base.pkgversion(ITensors.NDTensors)),
                     "KrylovKit"=>string(Base.pkgversion(KrylovKit))),
    "code"=>Dict("status"=>source_status, "before"=>SOURCE_BEFORE, "after"=>source_after),
    "acceptance_limits"=>limits, "points"=>points)
temporary, io = mktemp(dirname(output))
try
    TOML.print(io, record; sorted=true)
    close(io)
    mv(temporary, output; force=false)
finally
    isopen(io) && close(io)
    ispath(temporary) && rm(temporary)
end
println("Saved validation record: ", output)
record["status"] == "passed" || error("validation failed; all point outcomes were saved")
