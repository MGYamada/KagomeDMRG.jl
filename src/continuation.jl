"""
    FluxPolicy(; min_overlap, max_density_change, max_entropy_change,
                 max_schmidt_change, max_variance, max_truncation_error,
                 max_sweep_energy_change, max_cut_spread, consistency_tol=1e-9)

Explicit, finite diagnostic limits for one continuation experiment. The overlap
is an amplitude, not its square. Density, entropy and Schmidt-mean limits apply
between successive accepted points; cut spread uses cumulative right transfer
from the measured zero-flux state. These limits depend on geometry and accuracy
and are not universal convergence criteria or a topological branch certificate.
No expected pump value enters the policy.
"""
struct FluxPolicy
    min_overlap::Float64
    max_density_change::Float64
    max_entropy_change::Float64
    max_schmidt_change::Float64
    max_variance::Float64
    max_truncation_error::Float64
    max_sweep_energy_change::Float64
    max_cut_spread::Float64
    consistency_tol::Float64
    function FluxPolicy(; min_overlap, max_density_change, max_entropy_change,
                        max_schmidt_change, max_variance, max_truncation_error,
                        max_sweep_energy_change, max_cut_spread, consistency_tol=1e-9)
        values = Float64.((min_overlap, max_density_change, max_entropy_change,
            max_schmidt_change, max_variance, max_truncation_error,
            max_sweep_energy_change, max_cut_spread, consistency_tol))
        all(isfinite, values) && all(>=(0), values) && values[1] <= 1 &&
            values[end] > 0 || throw(ArgumentError("invalid finite flux policy limits"))
        new(values...)
    end
end

_flux_policy_record(policy) = Dict(string(k) => getfield(policy, k)
                                  for k in fieldnames(FluxPolicy))
_flux_data(point) = hasproperty(point, :diagnostics) ? point.diagnostics : point

function _flux_result(saved)
    return (; psi=saved.psi, sites=saved.sites, lattice=saved.lattice,
        hz=saved.hz, theta=saved.theta, gauge=saved.gauge, Q=saved.Q,
        settings=saved.settings, saved.diagnostics...)
end

function _continuation_solve(saved, lattice, theta; gauge, hz, outputlevel)
    s = saved.settings
    return run_dmrg(lattice, theta; sites=saved.sites, psi0=saved.psi,
        seed=s.seed, initial_linkdim=s.initial_linkdim, nsweeps=s.nsweeps,
        maxdim=s.maxdim, cutoff=s.cutoff, noise=s.noise, eigsolve_tol=s.eigsolve_tol,
        eigsolve_krylovdim=s.eigsolve_krylovdim, eigsolve_maxiter=s.eigsolve_maxiter,
        measure_variance=s.measure_variance, gauge, hz, outputlevel)
end

function _flux_schmidt(psi, bonds)
    return Dict(b => schmidt_diagnostics(psi, b) for b in bonds)
end

function _flux_diagnostics(point, previous, baseline, baseline_schmidt,
                           cuts, bonds, bulk_sites)
    data, prevdata = _flux_data(point), _flux_data(previous)
    schmidt = _flux_schmidt(point.psi, bonds)
    prevschmidt = _flux_schmidt(previous.psi, bonds)
    overlap = abs(inner(previous.psi, point.psi)) /
        sqrt(real(inner(previous.psi, previous.psi)) * real(inner(point.psi, point.psi)))
    transfers = Dict{String,Any}[]
    consistency = maximum(abs(schmidt[b].mean_left_sz - sum(data.sz[1:b])) for b in bonds)
    for transfer in spin_transfer(point.lattice, baseline.sz, data.sz; cuts)
        b = 3point.lattice.Ly * transfer.cut
        left = schmidt[b].mean_left_sz - baseline_schmidt[b].mean_left_sz
        consistency = max(consistency, abs(left - transfer.left), abs(-left - transfer.right))
        push!(transfers, Dict("cut"=>transfer.cut, "left"=>transfer.left,
            "right"=>transfer.right, "total"=>transfer.total,
            "schmidt_left"=>left, "schmidt_right"=>-left))
    end
    cut_spread = isempty(transfers) ? 0.0 :
        maximum(t["right"] for t in transfers) - minimum(t["right"] for t in transfers)
    # Upstream sweep energies are local Ritz values before the last truncation.
    # Include their difference from the measured final energy in this gate.
    sweep_change = max(abs(data.sweep_energies[end] - data.sweep_energies[end-1]),
                       abs(data.energy - data.sweep_energies[end]))
    return Dict{String,Any}(
        "overlap"=>Float64(overlap), "energy"=>data.energy,
        "variance"=>data.variance,
        "negative_variance_tolerance"=>100eps(Float64) * max(1.0, data.energy^2),
        "final_truncation_error"=>last(data.max_truncation_errors),
        "sweep_energy_change"=>sweep_change,
        "max_density_change"=>maximum(abs.(data.sz[bulk_sites] - prevdata.sz[bulk_sites])),
        "max_entropy_change"=>maximum(abs(schmidt[b].entropy-prevschmidt[b].entropy) for b in bonds),
        "max_schmidt_change"=>maximum(abs(schmidt[b].mean_left_sz-prevschmidt[b].mean_left_sz) for b in bonds),
        "cut_spread"=>cut_spread,
        "max_charge_error"=>max(abs(sum(data.sz)-point.Q/2),
            isempty(transfers) ? 0.0 : maximum(abs(t["total"]) for t in transfers)),
        "max_schmidt_consistency_error"=>consistency,
        "sz_profile"=>copy(data.sz), "sweep_energies"=>copy(data.sweep_energies),
        "measured_truncation_errors"=>copy(data.max_truncation_errors),
        "transfers"=>transfers,
        "schmidt"=>[Dict("bond"=>b, "probabilities"=>schmidt[b].probabilities,
            "left_q"=>schmidt[b].left_q, "entropy"=>schmidt[b].entropy,
            "mean_left_sz"=>schmidt[b].mean_left_sz,
            "variance_left_sz"=>schmidt[b].variance_left_sz) for b in bonds])
end

function _flux_reasons(diagnostics, policy)
    reasons = String[]
    keys = ("overlap", "energy", "variance", "negative_variance_tolerance",
        "final_truncation_error", "sweep_energy_change", "max_density_change",
        "max_entropy_change", "max_schmidt_change", "cut_spread",
        "max_charge_error", "max_schmidt_consistency_error")
    all(k -> diagnostics[k] isa Real && isfinite(diagnostics[k]), keys) ||
        return ["nonfinite_diagnostics"]
    diagnostics["overlap"] + policy.consistency_tol >= policy.min_overlap || push!(reasons, "low_overlap")
    diagnostics["overlap"] <= 1 + policy.consistency_tol || push!(reasons, "invalid_overlap")
    diagnostics["variance"] >= -diagnostics["negative_variance_tolerance"] ||
        push!(reasons, "negative_variance")
    abs(diagnostics["variance"]) <= policy.max_variance || push!(reasons, "variance")
    for (key, limit) in (("final_truncation_error", policy.max_truncation_error),
        ("sweep_energy_change", policy.max_sweep_energy_change),
        ("max_density_change", policy.max_density_change),
        ("max_entropy_change", policy.max_entropy_change),
        ("max_schmidt_change", policy.max_schmidt_change),
        ("cut_spread", policy.max_cut_spread),
        ("max_charge_error", policy.consistency_tol),
        ("max_schmidt_consistency_error", policy.consistency_tol))
        diagnostics[key] <= limit || push!(reasons, key)
    end
    return reasons
end

# Replace only this run's journal after closing the temporary file. Completed
# checkpoints remain immutable. This provides atomic visibility, not fsync.
function _write_flux_record(path, record)
    temporary, io = mktemp(dirname(path))
    try
        TOML.print(io, record; sorted=true)
        close(io)
        Base.Filesystem.rename(temporary, path)
    finally
        isopen(io) && close(io)
        ispath(temporary) && rm(temporary)
    end
    return nothing
end

"""
    continue_flux(lattice, targets; start, output_root, policy::FluxPolicy,
                  initial_step, min_step, max_trials=100, cuts=1:lattice.Lx-1,
                  diagnostic_bonds=..., bulk_sites=1:nsites(lattice),
                  gauge=:seam, hz=nothing, outputlevel=0)

Sequential, diagnostic-gated continuation through ordered, unwrapped targets.
`start` is a zero-flux `run_dmrg` result or an accepted checkpoint path. Every
trial reloads the last accepted checkpoint, rebuilds H, and retains its solver
settings and original zero-flux baseline. Rejections halve the step; accepted
steps retain that reduced size. Reversals are allowed. There is no automatic
increase of bond dimension, sweep count, or step size.

All policy limits are explicit. Variance measurement, at least two sweeps,
zero noise and seam gauge are required in this first driver. Raw overlaps in
the theta-dependent uniform gauge are not accepted as continuity diagnostics.
`bulk_sites` defaults to all sites and its actual selection is saved; users
must select physical bulk regions for larger cylinders. Geometric cuts are
always included among the Schmidt diagnostic bonds. No cuts exist at Lx=1.

Each run creates a unique directory under `output_root` with atomic TOML
journaling and separate trial/accepted snapshots. A diagnostic failure at the
minimum step, a numerical failure or the trial budget can return `:unresolved`.
If the initial diagnostics fail, no checkpoint is accepted and
`last_checkpoint` is `nothing`. `:completed` means these finite diagnostic
limits passed at the requested endpoints; it does not establish adiabaticity,
quantization, convergence across sizes, or phase identity. Restarting starts a
new journal with explicit new targets/policy/step and preserves the saved path.
"""
function continue_flux(lattice::KagomeCylinder, targets;
                       start, output_root::AbstractString, policy::FluxPolicy,
                       initial_step::Real, min_step::Real, max_trials::Integer=100,
                       cuts=collect(1:(lattice.Lx-1)),
                       diagnostic_bonds=unique([3lattice.Ly .* collect(cuts); div(nsites(lattice),2)]),
                       bulk_sites=collect(1:nsites(lattice)),
                       gauge::Symbol=:seam, hz=nothing, outputlevel::Integer=0,
                       _solver=_continuation_solve)
    gauge == :seam || throw(ArgumentError("continuation currently requires seam gauge"))
    endpoints = Float64.(collect(targets))
    !isempty(endpoints) && all(isfinite, endpoints) || throw(ArgumentError("targets must be nonempty and finite"))
    isfinite(initial_step) && isfinite(min_step) && 0 < min_step <= initial_step ||
        throw(ArgumentError("require finite 0 < min_step <= initial_step"))
    max_trials > 0 || throw(ArgumentError("max_trials must be positive"))
    N = nsites(lattice)
    c = Int.(collect(cuts))
    length(unique(c)) == length(c) && all(x -> 1 <= x < lattice.Lx, c) ||
        throw(ArgumentError("invalid geometric cuts"))
    lattice.Lx == 1 || !isempty(c) || throw(ArgumentError("select at least one geometric cut"))
    b = sort!(unique([Int.(collect(diagnostic_bonds)); 3lattice.Ly .* c]))
    !isempty(b) && all(x -> 1 <= x < N, b) || throw(ArgumentError("invalid diagnostic bonds"))
    bulk = Int.(collect(bulk_sites))
    !isempty(bulk) && length(unique(bulk)) == length(bulk) && all(x -> 1 <= x <= N, bulk) ||
        throw(ArgumentError("invalid monitored bulk sites"))
    config = _checkpoint_configuration(lattice, gauge, hz)
    fields = config["hz"]
    saved_start = start isa AbstractString ? load_checkpoint(start, lattice; gauge, hz=fields) : nothing
    initial = saved_start === nothing ? start : _flux_result(saved_start)
    saved_start === nothing && initial.theta != 0 &&
        throw(ArgumentError("a new continuation must start at zero flux"))
    _checkpoint_configuration(initial.lattice, initial.gauge, initial.hz) == config ||
        throw(ArgumentError("start model or gauge does not match continuation"))
    settings = _checkpoint_settings(initial.settings)
    settings.measure_variance && settings.nsweeps >= 2 && settings.noise == 0 ||
        throw(ArgumentError("continuation requires measured variance, at least two sweeps, and zero noise"))
    baseline = saved_start === nothing ? initial : saved_start.baseline
    theta_path = saved_start === nothing ? [0.0] : copy(saved_start.theta_path)
    mkpath(output_root)
    directory = mktempdir(abspath(output_root); prefix="trajectory-", cleanup=false)
    output_path = joinpath(directory, "trajectory.toml")
    accepted = String[]
    trials = Dict{String,Any}[]
    record = Dict{String,Any}("schema_version"=>1, "status"=>"running",
        "scope"=>"finite_diagnostic_gated_continuation", "phase_identification"=>"not_attempted",
        "configuration"=>config, "runtime"=>_checkpoint_runtime(),
        "provenance"=>_checkpoint_provenance(), "policy"=>_flux_policy_record(policy),
        "solver"=>Dict(string(k)=>v for (k,v) in pairs(settings)),
        "targets"=>endpoints, "initial_step"=>Float64(initial_step),
        "min_step"=>Float64(min_step), "max_trials"=>Int(max_trials),
        "cuts"=>c, "diagnostic_bonds"=>b, "bulk_sites"=>bulk,
        "axial_cut_coverage"=>!isempty(c), "restarted"=>saved_start !== nothing,
        "theta_path"=>theta_path, "accepted_checkpoints"=>String[], "trials"=>trials)
    last_checkpoint = nothing
    function finish(status, reason)
        record["status"] = string(status)
        record["reason"] = reason
        record["theta_path"] = copy(theta_path)
        record["last_checkpoint"] = last_checkpoint === nothing ? "" : relpath(last_checkpoint, directory)
        record["accepted_checkpoints"] = relpath.(accepted, Ref(directory))
        _write_flux_record(output_path, record)
        return (; status, reason, theta_path=copy(theta_path),
            accepted_checkpoints=copy(accepted), trials, last_checkpoint, output_path)
    end
    record["initial"] = Dict{String,Any}("theta"=>initial.theta, "status"=>"running",
        "checkpoint"=>"", "diagnostics"=>Dict{String,Any}(), "reasons"=>String[])
    _write_flux_record(output_path, record)
    # Saving validates state/model/energy/baseline before any acceptance label.
    # Invalid measured data must leave an unresolved record, not a stale run.
    initial_trial, initial_diag, baseline_schmidt = try
        snapshot = save_checkpoint(directory, initial; baseline, theta_path, status=:trial)
        record["initial"]["checkpoint"] = relpath(snapshot, directory)
        restored = load_checkpoint(snapshot, lattice; gauge, hz=fields, status=:trial)
        reference_schmidt = _flux_schmidt(restored.baseline.psi, b)
        diagnostic = _flux_diagnostics(restored, restored, restored.baseline,
                                       reference_schmidt, c, b, bulk)
        (snapshot, diagnostic, reference_schmidt)
    catch exception
        exception isa InterruptException && rethrow()
        exception isa OutOfMemoryError && rethrow()
        record["initial"]["status"] = "invalid"
        record["initial"]["reasons"] = ["initial_validation_error"]
        record["initial"]["exception_type"] = string(nameof(typeof(exception)))
        return finish(:unresolved, "initial_validation_error")
    end
    reasons = _flux_reasons(initial_diag, policy)
    record["initial"] = Dict("theta"=>initial.theta, "diagnostics"=>initial_diag,
        "checkpoint"=>relpath(initial_trial, directory), "reasons"=>reasons,
        "status"=>isempty(reasons) ? "accepted" : "rejected")
    isempty(reasons) || return finish(:unresolved, "initial_diagnostics")
    last_checkpoint = save_checkpoint(directory, initial; baseline, theta_path, status=:accepted)
    push!(accepted, last_checkpoint)
    record["accepted_checkpoints"] = relpath.(accepted, Ref(directory))
    record["last_checkpoint"] = relpath(last_checkpoint, directory)
    _write_flux_record(output_path, record)
    step = Float64(initial_step)
    for target in endpoints
        while last(theta_path) != target
            length(trials) < max_trials || return finish(:unresolved, "max_trials")
            from = last(theta_path)
            distance = abs(target-from)
            attempted_step = min(step, distance)
            theta = distance <= step ? target : from + copysign(attempted_step, target-from)
            theta != from || return finish(:unresolved, "step_underflow")
            saved = load_checkpoint(last_checkpoint, lattice; gauge, hz=fields,
                                    expected_theta=from)
            trial = Dict{String,Any}("theta"=>theta, "from_theta"=>from,
                "step"=>theta-from, "target"=>target, "status"=>"running",
                "checkpoint"=>"", "accepted_checkpoint"=>"",
                "reasons"=>String[], "diagnostics"=>Dict{String,Any}())
            push!(trials, trial)
            _write_flux_record(output_path, record)
            result = try
                # Keep the diagnostic reference independent of a mutating solver.
                _solver(deepcopy(saved), lattice, theta; gauge, hz=fields, outputlevel)
            catch exception
                exception isa InterruptException && rethrow()
                exception isa OutOfMemoryError && rethrow()
                trial["status"] = "solver_error"
                trial["reasons"] = ["solver_exception"]
                trial["exception_type"] = string(nameof(typeof(exception)))
                # Do not persist arbitrary exception strings or host details.
                if !(exception isa ErrorException || exception isa LinearAlgebra.LAPACKException)
                    return finish(:unresolved, "solver_exception")
                end
                nothing
            end
            if result !== nothing
                stage = "trial_validation_error"
                diagnostics = try
                    result.theta == theta || throw(ArgumentError("solver returned incorrect theta"))
                    _checkpoint_settings(result.settings) == merge(saved.settings, (initialization="provided_mps",)) ||
                        throw(ArgumentError("solver changed continuation settings"))
                    snapshot = save_checkpoint(directory, result; baseline=saved.baseline,
                        theta_path=[theta_path; theta], status=:trial)
                    trial["checkpoint"] = relpath(snapshot, directory)
                    stage = "diagnostic_error"
                    _flux_diagnostics(result, saved, saved.baseline,
                                      baseline_schmidt, c, b, bulk)
                catch exception
                    exception isa InterruptException && rethrow()
                    exception isa OutOfMemoryError && rethrow()
                    trial["status"] = "invalid_trial"
                    trial["reasons"] = [stage]
                    trial["exception_type"] = string(nameof(typeof(exception)))
                    return finish(:unresolved, stage)
                end
                trial["diagnostics"] = diagnostics
                trial["reasons"] = _flux_reasons(diagnostics, policy)
                if isempty(trial["reasons"])
                    last_checkpoint = save_checkpoint(directory, result; baseline=saved.baseline,
                        theta_path=[theta_path; theta], status=:accepted)
                    push!(theta_path, theta)
                    push!(accepted, last_checkpoint)
                    trial["status"] = "accepted"
                    trial["accepted_checkpoint"] = relpath(last_checkpoint, directory)
                    record["theta_path"] = copy(theta_path)
                    record["accepted_checkpoints"] = relpath.(accepted, Ref(directory))
                    record["last_checkpoint"] = relpath(last_checkpoint, directory)
                    _write_flux_record(output_path, record)
                    continue
                end
                trial["status"] = "rejected"
            end
            _write_flux_record(output_path, record)
            if attempted_step <= min_step
                return finish(:unresolved, "min_step")
            end
            step = max(Float64(min_step), attempted_step / 2)
        end
    end
    return finish(:completed, "targets_reached")
end
