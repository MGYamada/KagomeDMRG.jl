# These consistency limits and precision criteria come from the static research
# protocol. Measurement completion, consistency, and accuracy are separate axes.
const _STATIC_INTEGRITY_LIMITS = (
    norm=1e-12, charge=1e-12, zz_diagonal=1e-10, pm_diagonal=1e-10,
    zz_symmetry=1e-10, zz_imaginary=1e-10, pm_hermiticity=1e-10,
    fixed_charge_zz=1e-10, bond_sum=1e-9, schmidt_density=1e-10,
    schmidt_variance=1e-9, sz_remeasurement=1e-12)
const _STATIC_PRECISION_LIMITS = (
    energy_per_site=1e-6, sz_profile=1e-4, bond_profile=1e-4,
    variance_per_site=1e-5, truncation=1e-6)
const _STATIC_DIAGNOSTIC_FORMAT = "KagomeDMRG.static_diagnostics"

_static_dict(x) = Dict{String,Any}(string(k) => v for (k, v) in pairs(x))
_static_rows(x, f) = [f.(collect(row)) for row in eachrow(x)]

function _static_cuts(cuts, lattice)
    values = collect(cuts)
    all(c -> c isa Integer && !(c isa Bool) && 1 <= c < lattice.Lx, values) ||
        throw(ArgumentError("cuts must be integer cell boundaries in 1:Lx-1"))
    length(unique(values)) == length(values) || throw(ArgumentError("duplicate diagnostic cuts"))
    return Int.(values)
end

function _static_reference(reference, configuration, theta, identity)
    reference === nothing && return nothing
    reference isa AbstractDict || throw(ArgumentError("reference must be a static diagnostic report"))
    get(reference, "format", nothing) == _STATIC_DIAGNOSTIC_FORMAT &&
        get(reference, "schema_version", nothing) == 1 &&
        get(reference, "status", nothing) == "completed" &&
        get(get(reference, "integrity", Dict()), "status", nothing) == "passed" ||
        throw(ArgumentError("reference needs completed, consistent static diagnostics"))
    get(reference, "configuration", nothing) == configuration &&
        get(reference, "theta", nothing) == theta ||
        throw(ArgumentError("reference model, charge, gauge, or theta mismatch"))
    provenance = get(reference, "provenance", Dict())
    get(provenance, "source_sha256", nothing) == Dict(identity.source_sha256) &&
        get(provenance, "environment_sha256", nothing) == Dict(identity.environment_sha256) &&
        get(provenance, "manifest", nothing) == identity.manifest &&
        get(reference, "runtime", nothing) == Dict(identity.runtime) ||
        throw(ArgumentError("reference source, environment, or runtime mismatch"))
    measurements = reference["measurements"]
    for (key, count) in (("sz_profile", configuration["N"]),
                         ("bond_energy", length(configuration["bonds"])))
        values = measurements[key]
        values isa AbstractVector && length(values) == count &&
            all(v -> v isa Real && isfinite(v), values) ||
            throw(ArgumentError("invalid reference $key"))
    end
    measurements["energy"] isa Real && isfinite(measurements["energy"]) ||
        throw(ArgumentError("invalid reference energy"))
    _checkpoint_settings(reference["settings"])
    return reference
end

function _static_precision(measurements, settings, reference, integrity_passed, N)
    limits = _static_dict(_STATIC_PRECISION_LIMITS)
    values = Dict{String,Any}(
        "variance_per_site" => abs(measurements["variance_per_site"]),
        "truncation" => last(measurements["measured_truncation_errors"]))
    missing = String["energy_per_site", "sz_profile", "bond_profile"]
    comparison = "no_reference"
    if reference !== nothing
        previous = reference["measurements"]
        energies = measurements["sweep_energies"]
        changes = [abs(measurements["energy"] - previous["energy"]),
                   abs(measurements["energy"] - last(energies))]
        length(energies) >= 2 && push!(changes, abs(energies[end] - energies[end-1]))
        values["energy_per_site"] = maximum(changes)/N
        values["sz_profile"] = maximum(abs.(measurements["sz_profile"] - previous["sz_profile"]))
        values["bond_profile"] = maximum(abs.(measurements["bond_energy"] - previous["bond_energy"]))
        empty!(missing)
        previous_settings = _checkpoint_settings(reference["settings"])
        # A varying schedule is not a fixed-chi stationarity comparison, even if
        # two callers happened to supply identical varying schedules.
        fixed = all(==(last(settings.maxdim)), settings.maxdim) &&
                all(==(last(settings.maxdim)), previous_settings.maxdim)
        comparison = fixed ? "fixed_chi" : "chi_change"
    end
    passed = Dict(k => isfinite(v) && v <= limits[k] for (k, v) in values)
    status = !integrity_passed ? "not_evaluated" :
             !all(Base.values(passed)) ? "unmet" : isempty(missing) ? "passed" : "incomplete"
    return Dict{String,Any}("status" => status, "values" => values, "limits" => limits,
        "missing" => missing, "passed" => passed, "comparison" => comparison,
        "within_batch_sweep_change_measured" => length(measurements["sweep_energies"]) >= 2,
        "stationarity_evaluated" => integrity_passed && comparison == "fixed_chi",
        "stationarity_passed" => comparison == "fixed_chi" && status == "passed")
end

"""
    static_diagnostics(point; cuts=1:point.lattice.Lx-1, chirality=false,
                       reference=nothing, progress_callback=nothing)

Measure a completed `run_dmrg` or `load_checkpoint` point without optimizing it
or changing its MPS, settings, or stored diagnostics. Rebuild H from the actual
model, unwrapped flux, field, and gauge; recompute energy/density, correlations,
exchange bond energies, raw H†H variance, and Schmidt diagnostics at geometric
cell cuts. Optionally measure dressed triangle chirality. Empty cuts are valid.

Return a TOML-compatible schema-1 report with independent execution (`status`),
`integrity`, and `precision` statuses. The five static precision criteria keep
their existing thresholds. Without a prior same-model/theta `reference` report,
parent energy/density/bond changes are missing: precision is `incomplete` unless
a measured criterion is already `unmet`. A passed criterion set is not ground
state convergence or phase identification. A chi change is not stationarity.
One-sweep calls include available parent/solver changes without inventing a
within-batch sweep difference. Postprocessing never changes `measure_variance`.

`progress_callback(event)` receives scalar phase/status/elapsed notifications;
exceptions propagate. Source/environment identity is checked before and after.
"""
function static_diagnostics(point; cuts=1:point.lattice.Lx-1, chirality::Bool=false,
                            reference=nothing, progress_callback=nothing)
    started = time_ns()
    identity = _require_execution_identity(point)
    lattice = point.lattice
    N, Q = nsites(lattice), point.Q
    selected_cuts = _static_cuts(cuts, lattice)
    configuration = _checkpoint_configuration(lattice, point.gauge, point.hz; Q)
    theta = _finite_theta(point.theta)
    settings = _checkpoint_settings(point.settings)
    stored = _checkpoint_diagnostics(hasproperty(point, :diagnostics) ? point.diagnostics : point,
                                     settings, N, Q)
    reference = _static_reference(reference, configuration, theta, identity)
    reference_record = reference === nothing ? Dict{String,Any}("provided" => false) :
        Dict{String,Any}("provided" => true,
            "sha256" => bytes2hex(sha256(sprint(io -> TOML.print(io, reference; sorted=true)))))
    # Measurement routines may orthogonalize internally. Use a private copy even
    # when a backend currently promises to preserve its input.
    psi = deepcopy(point.psi)
    _checkpoint_state_check(psi, point.sites, Q)
    timings = Dict{String,Any}()
    function measure(f, name)
        _dmrg_progress(progress_callback, started; kind=:phase, phase=Symbol(name), status=:started)
        timings[name] = @elapsed value = f()
        _dmrg_progress(progress_callback, started; kind=:phase, phase=Symbol(name), status=:completed)
        return value
    end
    H = measure("hamiltonian") do
        twisted_exchange_mpo(point.sites, lattice, theta; gauge=point.gauge, hz=point.hz)
    end
    fresh = measure("state") do
        (; energy=inner(psi', H, psi), sz=sz_profile(psi), norm=norm(psi))
    end
    corr = measure("correlations") do
        spin_correlations(psi)
    end
    bonds = [b.Jz*real(corr.zz[b.i,b.j]) +
             b.Jxy*real(cis(bond_phase(lattice,b,theta; gauge=point.gauge))*corr.pm[b.i,b.j])
             for b in lattice.bonds]
    second_moment = measure("variance") do
        inner(H, psi, H, psi)
    end
    variance = real(second_moment) - stored.energy^2
    roundoff = 100eps(Float64)*max(1, stored.energy^2)
    schmidt = measure("schmidt") do
        map(selected_cuts) do cut
            b = 3lattice.Ly*cut
            d = schmidt_diagnostics(psi, b)
            merge(_static_dict(d), Dict("cut" => cut,
                "density_error" => abs(d.mean_left_sz-sum(stored.sz[1:b])),
                "variance_error" => abs(d.variance_left_sz -
                    (real(sum(corr.zz[1:b,1:b]))-sum(stored.sz[1:b])^2))))
        end
    end
    measurements = Dict{String,Any}(
        "energy" => stored.energy, "sz_profile" => copy(stored.sz), "bond_energy" => bonds,
        "column_sz" => [sum(stored.sz[3lattice.Ly*x+1:3lattice.Ly*(x+1)]) for x in 0:lattice.Lx-1],
        "correlations" => Dict("zz_real" => _static_rows(corr.zz, real),
            "zz_imag" => _static_rows(corr.zz, imag), "pm_real" => _static_rows(corr.pm, real),
            "pm_imag" => _static_rows(corr.pm, imag)),
        "HdaggerH_real" => real(second_moment), "HdaggerH_imag" => imag(second_moment),
        "variance" => variance, "variance_per_site" => variance/N,
        "variance_roundoff_scale" => roundoff, "schmidt" => schmidt,
        "remeasured_energy_real" => real(fresh.energy), "remeasured_energy_imag" => imag(fresh.energy),
        "sweep_energies" => copy(stored.sweep_energies),
        "measured_truncation_errors" => copy(stored.max_truncation_errors),
        "maxlinkdim" => maxlinkdim(psi), "chirality_measured" => chirality)
    if chirality
        measurements["chirality"] = measure("chirality") do
            triangle_chiralities(psi, lattice, theta; gauge=point.gauge)
        end
    end
    errors = Dict{String,Any}(
        "norm" => abs(fresh.norm-1), "charge" => abs(sum(fresh.sz)-Q/2),
        "energy_remeasurement" => abs(fresh.energy-stored.energy),
        "sz_remeasurement" => maximum(abs.(fresh.sz-stored.sz)),
        "zz_diagonal" => maximum(abs.(diag(corr.zz).-0.25)),
        "pm_diagonal" => maximum(abs.(diag(corr.pm).-(0.5 .+stored.sz))),
        "zz_symmetry" => maximum(abs.(corr.zz-transpose(corr.zz))),
        "zz_imaginary" => maximum(abs.(imag.(corr.zz))),
        "pm_hermiticity" => maximum(abs.(corr.pm-corr.pm')),
        "fixed_charge_zz" => maximum(abs.(vec(sum(corr.zz;dims=2))-(Q/2).*stored.sz)),
        "bond_sum" => abs(sum(bonds)-dot(point.hz,stored.sz)-stored.energy),
        "schmidt_density" => maximum((s["density_error"] for s in schmidt); init=0.0),
        "schmidt_variance" => maximum((s["variance_error"] for s in schmidt); init=0.0))
    limits = _static_dict(_STATIC_INTEGRITY_LIMITS)
    limits["energy_remeasurement"] = 1e-10*max(1,abs(stored.energy))
    failures = [k for (k,v) in errors if !isfinite(v) || v > limits[k]]
    isfinite(variance) && variance >= -roundoff || push!(failures,"invalid_variance")
    isfinite(second_moment) && abs(imag(second_moment)) <= roundoff ||
        push!(failures,"invalid_second_moment")
    all(isfinite,corr.zz) && all(isfinite,corr.pm) && all(isfinite,bonds) &&
        all(isfinite,fresh.sz) || push!(failures,"nonfinite_observables")
    all(v -> isfinite(v) && 0 <= v <= 1, stored.max_truncation_errors) ||
        push!(failures,"invalid_truncation")
    chirality && !all(isfinite, measurements["chirality"]) && push!(failures,"nonfinite_chirality")
    integrity_passed = isempty(failures)
    identity == _execution_identity() || error("implementation changed during diagnostics")
    return Dict{String,Any}("schema_version" => 1, "format" => _STATIC_DIAGNOSTIC_FORMAT,
        "status" => "completed", "configuration" => configuration, "theta" => theta, "Q" => Q,
        "settings" => _static_dict(settings), "cuts" => selected_cuts,
        "reference" => reference_record,
        "provenance" => _checkpoint_provenance(identity), "runtime" => Dict(identity.runtime),
        "measurements" => measurements, "measurement_seconds" => timings,
        "integrity" => Dict("status" => integrity_passed ? "passed" : "failed",
            "errors" => errors, "limits" => limits, "failures" => sort!(failures)),
        "precision" => _static_precision(measurements,settings,reference,integrity_passed,N))
end

_static_checkpoint_hashes(path) = Dict(name => _file_sha256(joinpath(path,name))
    for name in ("metadata.toml", "state.jls", "checksums.toml"))

"""
    diagnose_checkpoint(snapshot, lattice; output, status=:trial, ...)

Strictly load a trusted local checkpoint and call `static_diagnostics` into a
new `output` directory (its parent must exist). The model/load expectations
`Q`, `gauge`, `hz`, `expected_theta`, `expected_settings` have the same meaning
as in `load_checkpoint`. Measurement options are `cuts`, `chirality`, `reference`.

Write `diagnostics.toml` atomically, first `running`, then `completed` or `failed`.
Load/measurement exceptions are rethrown after a failed record is written; only
the exception type and phase are saved, never arbitrary exception text. A killed
process may leave `running`, which is not success. Reruns require a fresh output;
existing directories and outputs inside the snapshot are rejected. All three
checkpoint hashes are checked before/after, including on failure. Parent files
and solver settings are never rewritten. Completion alone says nothing about
integrity or precision; inspect both fields of the returned report.
"""
function diagnose_checkpoint(snapshot::AbstractString, lattice::KagomeCylinder;
        output::AbstractString, Q=nothing, gauge::Symbol=:seam, hz=nothing,
        expected_theta=nothing, expected_settings=nothing, status::Symbol=:trial,
        cuts=1:lattice.Lx-1, chirality::Bool=false, reference=nothing)
    parent = realpath(snapshot)
    destination = joinpath(realpath(dirname(abspath(output))), basename(abspath(output)))
    relative = relpath(destination,parent)
    first(splitpath(relative)) == ".." ||
        throw(ArgumentError("diagnostic output must be outside the checkpoint"))
    ispath(destination) && throw(ArgumentError("diagnostic output already exists"))
    mkdir(destination)
    report = Dict{String,Any}("schema_version" => 1, "format" => _STATIC_DIAGNOSTIC_FORMAT,
        "status" => "running", "active_phase" => "checkpoint_load",
        "integrity" => Dict("status" => "not_evaluated"),
        "precision" => Dict("status" => "not_evaluated"))
    path = joinpath(destination,"diagnostics.toml")
    persist() = _write_flux_record(path,report)
    persist()
    hashes = nothing
    try
        hashes = _static_checkpoint_hashes(parent)
        report["checkpoint"] = Dict("path" => relpath(parent,destination), "sha256" => hashes)
        persist()
        saved = load_checkpoint(parent,lattice; Q,gauge,hz,expected_theta,expected_settings,status)
        function progress(event)
            report["active_phase"] = string(event.phase)
            report["phase_status"] = string(event.status)
            persist()
        end
        measured = static_diagnostics(saved; cuts,chirality,reference,progress_callback=progress)
        merge!(report,measured)
        report["active_phase"] = "completed"
        report["phase_status"] = "completed"
    catch exception
        report["status"] = "failed"
        report["exception_type"] = string(nameof(typeof(exception)))
        rethrow()
    finally
        unchanged = hashes !== nothing && try
            _static_checkpoint_hashes(parent) == hashes
        catch
            false
        end
        report["checkpoint_unchanged"] = unchanged
        if !unchanged
            report["status"] = "failed"
            report["input_status"] = "changed_or_unreadable"
            report["integrity"]["status"] = "not_evaluated"
            report["precision"]["status"] = "not_evaluated"
        end
        persist()
    end
    report["checkpoint_unchanged"] || error("checkpoint changed during diagnostics")
    return report
end
