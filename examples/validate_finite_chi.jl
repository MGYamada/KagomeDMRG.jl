#!/usr/bin/env julia
# Calibrate the patched factorization and compare finite-chi single-point
# optimizations. These points are not accepted flux-continuation trajectories.
# julia --project=. --startup-file=no --threads=1 examples/validate_finite_chi.jl NEW_OUTPUT [UPSTREAM_CAPTURE [CHI:SWEEPS,...]]
using KagomeDMRG, ITensors, ITensorMPS, LinearAlgebra, Serialization, SHA, TOML, Dates
BLAS.set_num_threads(1)
ITensors.disable_threaded_blocksparse()
const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "test", "reference_ed.jl"))
include(joinpath(ROOT, "test", "reference_eigensolve.jl"))
include(joinpath(ROOT, "test", "itensor_helpers.jl"))

length(ARGS) in 1:3 || error("usage: validate_finite_chi.jl new-output-directory [trusted-upstream-capture [chi:sweeps,...]]")
const ALL_CASES = ((128,6), (128,12), (256,6), (256,12), (512,6))
const CASES = length(ARGS) == 3 ? Tuple(Tuple(parse.(Int, split(case, ':')))
    for case in split(ARGS[3], ',')) : ALL_CASES
!isempty(CASES) && all(case in ALL_CASES for case in CASES) &&
    length(unique(CASES)) == length(CASES) || error("select unique cases from the fixed calibration protocol")
const OUTPUT = abspath(ARGS[1])
ispath(OUTPUT) && error("choose a new output directory")
mkpath(OUTPUT)
const CODE = KagomeDMRG._checkpoint_provenance()
const EXTRA_SOURCES = ("examples/validate_finite_chi.jl", "test/reference_ed.jl",
    "test/reference_eigensolve.jl", "test/itensor_helpers.jl")
extra_hashes() = Dict(p => bytes2hex(sha256(read(joinpath(ROOT, p)))) for p in EXTRA_SOURCES)
const EXTRA_BEFORE = extra_hashes()
const SOLVER = (; cutoff=0.0, noise=0.0, eigsolve_tol=1e-11,
    eigsolve_krylovdim=40, eigsolve_maxiter=20, measure_variance=true)
const ED_LIMITS = Dict("energy_error_per_site"=>1e-8, "residual_norm"=>2e-6,
    "max_sz_error"=>1e-6, "ground_space_leakage"=>2e-5)
const RECORD = Dict{String,Any}("schema_version"=>1, "status"=>"running",
    "recorded_at_utc"=>string(now(UTC)),
    "scope"=>"finite_chi_accuracy_and_truncation_calibration",
    "phase_identification"=>"not_attempted", "quantization"=>"not_tested",
    "continuation_acceptance"=>"not_attempted_single_point_optimizations",
    "code"=>CODE, "extra_source_sha256"=>EXTRA_BEFORE,
    "runtime"=>merge(KagomeDMRG._checkpoint_runtime(), Dict(
        "julia_threads"=>Threads.nthreads(), "blas_threads"=>BLAS.get_num_threads())),
    "ed_limits"=>ED_LIMITS, "solver"=>Dict(string(k)=>v for (k,v) in pairs(SOLVER)),
    "planned_cases"=>[Dict("chi"=>chi, "nsweeps"=>sweeps) for (chi,sweeps) in CASES],
    "references"=>Dict{String,Any}[],
    "runs"=>Dict{String,Any}[], "boundary_replays"=>Dict{String,Any}[])
persist() = KagomeDMRG._write_flux_record(joinpath(OUTPUT, "validation.toml"), RECORD)

function projection_loss(input, kept)
    coefficient = dot(vec(kept), vec(input)) / sum(abs2, kept)
    return sum(abs2, input - coefficient * kept) / sum(abs2, input)
end

# Inspect only the central bond, on both sweep directions. This tracing method
# in the study process delegates the update to unchanged ITensorMPS code.
const AUDIT = Ref{Any}(nothing)
const UPDATE = Ref(0)
function ITensorMPS.replacebond!(PH::ITensorMPS.ProjMPO, psi::MPS,
                                bond::Int, phi::ITensor; kwargs...)
    UPDATE[] += 1
    monitored = AUDIT[] !== nothing && bond == 9
    left = monitored ? collect(commoninds(phi, psi[bond])) : nothing
    right = monitored ? collect(uniqueinds(phi, left)) : nothing
    input = monitored ? reshape(Array(phi, left..., right...), prod(dim, left), :) : nothing
    spec = ITensorMPS.replacebond!(psi, bond, phi; kwargs...)
    if monitored
        kept = reshape(Array(psi[bond] * psi[bond+1], left..., right...), size(input))
        rank = dim(linkind(psi, bond))
        weights = svdvals(input).^2
        actual_loss = projection_loss(input, kept)
        best_loss = sum(weights[(rank+1):end]) / sum(weights)
        row = Dict{String,Any}("update"=>UPDATE[], "bond"=>bond,
            "sweep"=>div(UPDATE[]-1, 2(length(psi)-1))+1, "ortho"=>kwargs[:ortho],
            "retained_rank"=>rank, "reported_rank"=>length(eigs(spec)),
            "input_norm_squared"=>sum(weights), "reported_loss"=>truncerror(spec),
            "dense_best_rank_loss"=>best_loss, "measured_projection_loss"=>actual_loss,
            "reported_loss_error"=>abs(truncerror(spec)-actual_loss),
            "best_rank_loss_error"=>abs(best_loss-actual_loss))
        row["passed"] = rank == length(eigs(spec)) && rank <= kwargs[:maxdim] &&
            row["reported_loss_error"] < 1e-11 && row["best_rank_loss_error"] < 1e-11
        push!(AUDIT[], row)
    end
    return spec
end

function replay_boundary(directory)
    metadata = TOML.parsefile(joinpath(directory, "capture.toml"))
    metadata["status"] == "captured_guard_failure" && metadata["sources_unchanged"] ||
        error("unverified boundary capture")
    file = joinpath(directory, "factorization_input.jls")
    bytes2hex(sha256(read(file))) == metadata["payload_sha256"] || error("boundary checksum mismatch")
    # Explicitly authorized replay of the trusted numerical input across a
    # function-only backend patch. This is not a production checkpoint restart.
    captured = open(deserialize, file)
    captured.kwargs.maxdim == 256 && captured.kwargs.cutoff == 0.0 &&
        get(captured.kwargs, :eigen_perturbation, nothing) === nothing &&
        metadata["solver"]["noise"] == 0.0 ||
        error("boundary replay requires the calibrated rank-256 noiseless capture")
    phi, left = captured.phi, captured.left
    ITensors.checkflux(phi)
    val(flux(phi), "Sz") == metadata["input_charge"] || error("boundary charge mismatch")
    right = collect(uniqueinds(phi, left))
    input = reshape(Array(phi, left..., right...), prod(dim, left), :)
    weights = svdvals(input).^2
    reference_loss = sum(weights[257:end]) / sum(weights)
    for decomposition in ("svd", "eigen"), ortho in ("left", "right")
        L, R, spec = factorize(phi, left; which_decomp=decomposition,
            ortho, maxdim=256, cutoff=0.0)
        kept = reshape(Array(L * R, left..., right...), size(input))
        loss = projection_loss(input, kept)
        rank = dim(commonind(L, R))
        row = Dict{String,Any}("decomposition"=>decomposition, "ortho"=>ortho,
            "retained_rank"=>rank, "reported_rank"=>length(eigs(spec)),
            "measured_projection_loss"=>loss, "dense_best_rank_loss"=>reference_loss,
            "reported_loss"=>truncerror(spec),
            "input_payload_sha256"=>metadata["payload_sha256"],
            "upstream_retained_rank"=>metadata["actual_rank"],
            "upstream_reported_rank"=>metadata["reported_rank"],
            "upstream_reported_loss"=>metadata["reported_truncation_error"])
        row["passed"] = rank == length(eigs(spec)) == 256 &&
            abs(loss-reference_loss) < 1e-12 && abs(loss-truncerror(spec)) < 1e-12
        push!(RECORD["boundary_replays"], row)
    end
    persist()
    all(r["passed"] for r in RECORD["boundary_replays"]) || error("boundary replay failed")
end

const REFERENCES = Dict{Float64,Any}()
function reference(theta)
    return get!(REFERENCES, theta) do
        elapsed = @elapsed ed = reference_eigensystem(2, 3, theta; seed=5678, blocksize=4)
        info = Dict{String,Any}("theta"=>theta, "elapsed_seconds"=>elapsed,
            "values"=>ed.values, "residual_norms"=>ed.residual_norms,
            "orthogonality_error"=>ed.orthogonality_error, "converged"=>ed.converged,
            "settings"=>Dict(string(k)=>v for (k,v) in pairs(ed.settings)),
            "solver_info"=>Dict(string(k)=>v for (k,v) in pairs(ed.solver_info)))
        push!(RECORD["references"], info)
        persist()
        ed.converged || error("independent reference did not converge")
        last(ed.values)-first(ed.values) > 1e-8 || error("ground cluster unresolved")
        open(joinpath(OUTPUT, "reference_$(length(REFERENCES)+1).jls"), "w") do io
            serialize(io, (; theta, basis=ed.basis, values=ed.values, vectors=ed.vectors))
        end
        return ed
    end
end

function metrics(result, baseline)
    ed = reference(result.theta)
    v = sector_amplitudes(result.psi, result.sites, ed.basis)
    cluster = findall(x -> x-first(ed.values) < 1e-8, ed.values)
    projected = ed.vectors[:,cluster] * (ed.vectors[:,cluster]' * v)
    residual = norm(ed.H*v-result.energy*v)
    schmidt = schmidt_diagnostics(result.psi, 9)
    base_schmidt = schmidt_diagnostics(baseline.psi, 9)
    transfer = only(spin_transfer(result.lattice, baseline.sz, result.sz))
    row = Dict{String,Any}("theta"=>result.theta, "energy"=>result.energy,
        "energy_error_per_site"=>abs(result.energy-first(ed.values))/18,
        "residual_norm"=>residual, "raw_variance"=>result.variance,
        "variance_residual_squared_error"=>abs(result.variance-residual^2),
        "ground_space_leakage"=>norm(v-projected),
        "max_sz_error"=>maximum(abs.(result.sz-reference_sz(projected,ed.basis,18))),
        "sz_profile"=>result.sz, "total_sz_error"=>abs(sum(result.sz)-1),
        "state_norm_error"=>abs(norm(v)-1), "maxlinkdim"=>maxlinkdim(result.psi),
        "sweep_energies"=>result.sweep_energies,
        "measured_truncation_errors"=>result.max_truncation_errors,
        "schmidt_entropy"=>schmidt.entropy, "schmidt_mean_left_sz"=>schmidt.mean_left_sz,
        "schmidt_density_error"=>abs(schmidt.mean_left_sz-sum(result.sz[1:9])),
        "transfer_right"=>transfer.right, "transfer_left"=>transfer.left,
        "transfer_total"=>transfer.total,
        "schmidt_transfer_right"=>-(schmidt.mean_left_sz-base_schmidt.mean_left_sz))
    row["ed_comparison"] = all(row[k] <= limit for (k,limit) in ED_LIMITS) ? "within_limits" : "outside_limits"
    return row
end

persist()
try
    length(ARGS) >= 2 && replay_boundary(abspath(ARGS[2]))
    lattice = kagome_cylinder(2, 3)
    RECORD["configuration"] = KagomeDMRG._checkpoint_configuration(lattice, :seam, nothing)
    # A common reference failure invalidates the entire comparison, so resolve
    # it before spending time on any of the finite-chi optimizations.
    reference(0.0)
    reference(0.37)
    # These cases are fixed before inspecting outcomes. Flux points in each
    # case remain sequential, using its own zero-flux state and measured profile.
    for (chi, nsweeps) in CASES
        run = Dict{String,Any}("chi"=>chi, "nsweeps"=>nsweeps, "seed"=>11,
            "status"=>"running", "points"=>Dict{String,Any}[])
        push!(RECORD["runs"], run)
        persist()
        baseline = nothing
        try
            for theta in (0.0, 0.37)
                audit = Dict{String,Any}[]
                AUDIT[] = audit
                UPDATE[] = 0
                elapsed = @elapsed result = if baseline === nothing
                    run_dmrg(lattice, theta; seed=11, nsweeps, maxdim=chi, SOLVER...)
                else
                    run_dmrg(lattice, theta; sites=baseline.sites, psi0=baseline.psi,
                        seed=11, nsweeps, maxdim=chi, SOLVER...)
                end
                AUDIT[] = nothing
                baseline === nothing && (baseline = result)
                snapshot = save_checkpoint(OUTPUT, result; baseline,
                    theta_path=theta == 0 ? [0.0] : [0.0, theta], status=:trial)
                row = metrics(result, baseline)
                row["checkpoint"] = relpath(snapshot, OUTPUT)
                row["elapsed_seconds"] = elapsed
                row["central_update_calibration"] = audit
                row["calibration_passed"] = length(audit) == 2nsweeps && all(x["passed"] for x in audit)
                push!(run["points"], row)
                persist()
                println("chi=", chi, " sweeps=", nsweeps, " theta=", theta,
                    " energy_error/site=", row["energy_error_per_site"],
                    " residual=", row["residual_norm"],
                    " calibration=", row["calibration_passed"])
            end
            run["status"] = "completed_with_recorded_accuracy"
        catch exception
            AUDIT[] = nothing
            exception isa InterruptException && rethrow()
            exception isa OutOfMemoryError && rethrow()
            run["status"] = "error"
            run["exception_type"] = string(nameof(typeof(exception)))
            showerror(stderr, exception)
            println(stderr)
        end
        persist()
    end
    complete = all(r["status"] == "completed_with_recorded_accuracy" &&
        all(p["calibration_passed"] for p in r["points"]) for r in RECORD["runs"])
    RECORD["status"] = complete ? "completed_with_recorded_accuracy" : "incomplete_or_calibration_failure"
catch exception
    exception isa InterruptException && rethrow()
    exception isa OutOfMemoryError && rethrow()
    RECORD["status"] = "error"
    RECORD["exception_type"] = string(nameof(typeof(exception)))
    showerror(stderr, exception)
    println(stderr)
end
RECORD["sources_unchanged"] = try
    extra_hashes() == EXTRA_BEFORE &&
        KagomeDMRG._checkpoint_provenance()["source_sha256"] == CODE["source_sha256"]
catch exception
    exception isa InterruptException && rethrow()
    exception isa OutOfMemoryError && rethrow()
    false
end
RECORD["sources_unchanged"] || (RECORD["status"] = "source_changed")
persist()
println("Saved finite-chi validation: ", joinpath(OUTPUT, "validation.toml"))
RECORD["status"] == "completed_with_recorded_accuracy" || error("study did not complete; outcomes retained")
