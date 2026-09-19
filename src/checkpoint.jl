# Schema 1 is a trusted, local restart format, not a long-term interchange
# format. Check primitive TOML metadata and checksums before deserializing.
const _CHECKPOINT_FORMAT = "KagomeDMRG.local_checkpoint"

_checkpoint_require(condition, message) = condition || throw(ArgumentError(message))

function _checkpoint_charge(N, stored_Q; Q=nothing)
    _checkpoint_require(stored_Q isa Integer && !(stored_Q isa Bool),
                        "checkpoint charge must be an integer")
    charge = target_sector(N; Q=stored_Q).Q
    Q === nothing || _checkpoint_require(target_sector(N; Q).Q == charge,
                                         "requested charge does not match checkpoint")
    return charge
end

function _checkpoint_configuration(lattice::KagomeCylinder, gauge, hz; Q=nothing)
    sector = target_sector(nsites(lattice); Q)
    _checkpoint_require(lattice.Lx >= 1 && lattice.Ly >= 3 &&
        sector.N == 3lattice.Lx*lattice.Ly, "invalid checkpoint geometry")
    _checkpoint_require(gauge in (:seam, :uniform), "unknown checkpoint gauge")
    fields = hz === nothing ? zeros(sector.N) : Float64.(hz)
    _checkpoint_require(length(fields) == sector.N && all(isfinite, fields),
                        "invalid checkpoint longitudinal fields")
    _checkpoint_require(all(s.index == i && all(isfinite, s.position)
        for (i, s) in enumerate(lattice.sites)), "invalid checkpoint site ordering")
    canonical = kagome_cylinder(lattice.Lx, lattice.Ly)
    _checkpoint_require(lattice.sites == canonical.sites,
                        "checkpoint sites must match the documented kagome geometry")
    _checkpoint_require(all(1 <= b.i <= sector.N && 1 <= b.j <= sector.N &&
        b.i != b.j && isfinite(b.Jxy) && isfinite(b.Jz) for b in lattice.bonds),
        "invalid checkpoint bonds")
    families = bond_families(lattice)
    return Dict{String,Any}(
        "Lx" => lattice.Lx, "Ly" => lattice.Ly, "N" => sector.N, "Q" => sector.Q,
        "charge_convention" => "q=2Sz", "ordering" => "x_then_y_then_A_B_C",
        "axis_boundary" => "open", "circumference_boundary" => "periodic",
        "termination" => "complete_cells_with_outgoing_open_axis_bonds_removed",
        "wrap" => [lattice.Ly/2, lattice.Ly*sqrt(3)/2],
        "exchange_phase" => "S+_i_S-_j_has_exp(+im*bond_phase)",
        "bond_family_convention" => "J1_nn_J2_second_J3_hexagon_opposite_other_unclassified",
        "gauge" => string(gauge), "hz" => fields,
        "sites" => [Dict("index" => s.index, "x" => s.x, "y" => s.y,
            "sublattice" => string(s.sublattice), "position" => collect(s.position))
            for s in lattice.sites],
        "bonds" => [Dict("i" => b.i, "j" => b.j, "Jxy" => b.Jxy,
            "Jz" => b.Jz, "wy" => b.wy, "family" => string(families[k]))
            for (k, b) in enumerate(lattice.bonds)])
end

function _checkpoint_site_record(s)
    return Dict("id" => string(id(s)), "prime_level" => plev(s),
        "direction" => Int(dir(s)), "tags" => string(tags(s)), "dim" => dim(s),
        "q" => [val(qn(s => k), "Sz") for k in 1:dim(s)])
end

function _checkpoint_state_check(psi, sites, Q)
    _checkpoint_require(psi isa MPS, "checkpoint state must be an MPS")
    _check_sites(sites, length(psi))
    _checkpoint_require(length(unique(id.(sites))) == length(sites),
                        "checkpoint site IDs must be unique")
    _checkpoint_require(all(s -> dir(s) == ITensors.Out && plev(s) == 0, sites),
                        "checkpoint physical indices must be outgoing and unprimed")
    for i in eachindex(sites)
        _checkpoint_require(_checkpoint_site_record(siteind(psi, i)) ==
            _checkpoint_site_record(sites[i]), "checkpoint MPS site indices do not match")
        _checkpoint_require(eltype(psi[i]) == ComplexF64 && isfinite(norm(psi[i])),
                            "checkpoint MPS must have finite ComplexF64 tensors")
        _checkpoint_require(order(psi[i]) == (i in (1, length(psi)) ? 2 : 3),
                            "checkpoint must be an open-boundary MPS")
        ITensors.checkflux(psi[i])
    end
    _checkpoint_require(all(i -> length(commoninds(psi[i], psi[i+1])) == 1,
                            1:(length(psi)-1)), "checkpoint MPS links do not match")
    _checkpoint_require(flux(psi) == QN("Sz", Q), "checkpoint total charge mismatch")
    norm_squared = inner(psi, psi)
    _checkpoint_require(isfinite(norm_squared) && abs(norm_squared - 1) <= 1e-10,
                        "checkpoint state must be normalized")
    return nothing
end

_checkpoint_value(settings, key) = settings isa AbstractDict ?
    settings[string(key)] : getproperty(settings, key)

function _checkpoint_settings(settings)
    get(key) = _checkpoint_value(settings, key)
    s = (; seed=Int(get(:seed)), initial_linkdim=Int(get(:initial_linkdim)),
        nsweeps=Int(get(:nsweeps)), maxdim=Int.(get(:maxdim)),
        cutoff=Float64(get(:cutoff)), noise=Float64(get(:noise)),
        eigsolve_tol=Float64(get(:eigsolve_tol)),
        eigsolve_krylovdim=Int(get(:eigsolve_krylovdim)),
        eigsolve_maxiter=Int(get(:eigsolve_maxiter)),
        measure_variance=Bool(get(:measure_variance)),
        initialization=String(get(:initialization)))
    _checkpoint_require(s.initial_linkdim > 0 && s.nsweeps > 0 &&
        !isempty(s.maxdim) && all(>(0), s.maxdim) &&
        isfinite(s.cutoff) && s.cutoff >= 0 && isfinite(s.noise) && s.noise >= 0 &&
        isfinite(s.eigsolve_tol) && s.eigsolve_tol > 0 && s.eigsolve_krylovdim >= 2 &&
        s.eigsolve_maxiter > 0 && s.initialization in ("random_fixed_charge", "provided_mps"),
        "invalid checkpoint solver settings")
    return s
end

function _checkpoint_profile(sz, N, Q)
    _checkpoint_require(length(sz) == N && all(isfinite, sz) &&
        all(x -> abs(x) <= 0.5 + 1e-10, sz) && abs(sum(sz) - Q/2) <= 1e-9,
        "invalid checkpoint spin profile")
    return Float64.(sz)
end

function _checkpoint_diagnostics(state, settings, N, Q)
    get(key) = _checkpoint_value(state, key)
    sz = _checkpoint_profile(get(:sz), N, Q)
    energy, local_energy = Float64(get(:energy)), Float64(get(:local_energy))
    _checkpoint_require(isfinite(energy) && isfinite(local_energy),
                        "checkpoint energies must be finite")
    energies = Float64.(get(:sweep_energies))
    errors = Float64.(get(:max_truncation_errors))
    _checkpoint_require(length(energies) == length(errors) == settings.nsweeps &&
        all(isfinite, energies) && all(x -> isfinite(x) && x >= 0, errors),
        "invalid checkpoint sweep diagnostics")
    measured = state isa AbstractDict ? state["variance_measured"] : get(:variance) !== nothing
    _checkpoint_require(measured == settings.measure_variance,
                        "checkpoint variance measurement setting mismatch")
    variance = measured ? Float64(get(:variance)) : nothing
    _checkpoint_require(variance === nothing || isfinite(variance),
                        "checkpoint variance must be finite when measured")
    return (; energy, local_energy, variance, sz, sweep_energies=energies,
            max_truncation_errors=errors)
end

function _checkpoint_flux_path(theta_path, theta)
    values = Float64.(collect(theta_path))
    _checkpoint_require(!isempty(values) && all(isfinite, values) &&
        first(values) == 0 && last(values) == theta, "checkpoint flux path must start at zero and end at theta")
    return values
end

function _checkpoint_energy_check(psi, sites, lattice, theta, gauge, hz, energy)
    H = twisted_exchange_mpo(sites, lattice, theta; gauge, hz)
    measured = inner(psi', H, psi)
    _checkpoint_require(isfinite(measured) &&
        abs(measured - energy) <= 1e-10 * max(1, abs(energy)),
        "checkpoint energy does not match its MPS and Hamiltonian")
    return nothing
end

"""
    save_checkpoint(root, result; baseline, theta_path, status=:trial)

Save an immutable, completed-point MPS snapshot below `root/trial/` or
`root/accepted/` and return its directory. `baseline` is a measured zero-flux
`run_dmrg` result with the same model, total charge and site indices, or the `baseline`
returned by `load_checkpoint`. Its zero-flux MPS is also saved so the measured
baseline can be checked again on load. `theta_path` starts at zero and ends at the
unwrapped `result.theta`; it may reverse direction. The caller supplies the
path and acceptance decision; this function does not certify branch continuity.

The directory becomes visible only after all files have been closed and a
same-parent rename succeeds. Existing snapshots are never replaced. This
guarantees atomic visibility, not power-loss durability. The Julia Serialization
payload is for trusted local use with the same runtime, source, and dependencies;
it is not a portable or long-term archival format. No sweep environments are
saved. Restart begins a new DMRG batch from the completed MPS.
"""
function save_checkpoint(root::AbstractString, result; baseline, theta_path,
                         status::Symbol=:trial)
    _checkpoint_require(status in (:trial, :accepted), "unknown checkpoint status")
    execution_identity = _require_execution_identity(result)
    _require_execution_identity(baseline, execution_identity)
    charge = _checkpoint_charge(nsites(result.lattice), result.Q)
    config = _checkpoint_configuration(result.lattice, result.gauge, result.hz; Q=charge)
    N, Q = config["N"], config["Q"]
    _checkpoint_state_check(result.psi, result.sites, Q)
    settings = _checkpoint_settings(result.settings)
    diagnostics = _checkpoint_diagnostics(result, settings, N, Q)
    _checkpoint_require(isapprox(sz_profile(result.psi), diagnostics.sz; atol=1e-10, rtol=0),
                        "checkpoint profile does not match the MPS")
    theta = _finite_theta(result.theta)
    _checkpoint_energy_check(result.psi, result.sites, result.lattice, theta,
                             result.gauge, result.hz, diagnostics.energy)
    path = _checkpoint_flux_path(theta_path, theta)
    _checkpoint_require(hasproperty(baseline, :psi),
                        "baseline must include its measured zero-flux MPS")
    baseline_Q = _checkpoint_charge(nsites(baseline.lattice), baseline.Q)
    _checkpoint_require(baseline.theta == 0 && baseline_Q == Q &&
        _checkpoint_configuration(baseline.lattice, baseline.gauge, baseline.hz; Q=baseline_Q) == config &&
        _checkpoint_site_record.(baseline.sites) == _checkpoint_site_record.(result.sites),
        "zero-flux baseline configuration or site indices do not match")
    baseline_sz = _checkpoint_profile(baseline.sz, N, Q)
    _checkpoint_require(length(path) != 1 ||
        isapprox(diagnostics.sz, baseline_sz; atol=1e-10, rtol=0),
        "initial checkpoint must use its own measured zero-flux baseline")
    _checkpoint_state_check(baseline.psi, baseline.sites, Q)
    _checkpoint_require(isapprox(sz_profile(baseline.psi), baseline_sz; atol=1e-10, rtol=0),
                        "baseline profile does not match the zero-flux MPS")
    state = Dict{String,Any}("theta" => theta, "Q" => Q,
        "site_indices" => _checkpoint_site_record.(result.sites),
        "sz" => diagnostics.sz, "energy" => diagnostics.energy,
        "local_energy" => diagnostics.local_energy,
        "variance_measured" => diagnostics.variance !== nothing,
        "sweep_energies" => diagnostics.sweep_energies,
        "max_truncation_errors" => diagnostics.max_truncation_errors)
    diagnostics.variance === nothing || (state["variance"] = diagnostics.variance)
    provenance = _checkpoint_provenance(execution_identity)
    metadata = Dict{String,Any}("schema_version" => 1, "format" => _CHECKPOINT_FORMAT,
        "status" => string(status), "created_unix" => time(),
        "runtime" => Dict(execution_identity.runtime), "provenance" => provenance,
        "configuration" => config, "state" => state,
        "settings" => Dict(string(k) => v for (k,v) in pairs(settings)),
        "baseline" => Dict("theta" => 0.0, "sz" => baseline_sz), "theta_path" => path)
    parent = joinpath(abspath(root), string(status))
    mkpath(parent)
    temporary = mktempdir(parent; prefix=".pending-", cleanup=false)
    destination = joinpath(parent, replace(basename(temporary), ".pending-" => "checkpoint-"))
    try
        open(joinpath(temporary, "metadata.toml"), "w") do io
            TOML.print(io, metadata; sorted=true)
        end
        open(joinpath(temporary, "state.jls"), "w") do io
            serialize(io, (; psi=result.psi, sites=result.sites, baseline_psi=baseline.psi))
        end
        checksums = Dict("schema_version" => 1, "files" => Dict(
            name => Dict("sha256" => _file_sha256(joinpath(temporary, name)),
                         "bytes" => filesize(joinpath(temporary, name)))
            for name in ("metadata.toml", "state.jls")))
        open(joinpath(temporary, "checksums.toml"), "w") do io
            TOML.print(io, checksums; sorted=true)
        end
        _checkpoint_require(execution_identity == _execution_identity(),
                            "source changed while saving checkpoint")
        _checkpoint_require(!ispath(destination), "checkpoint destination already exists")
        Base.Filesystem.rename(temporary, destination)
    finally
        isdir(temporary) && rm(temporary; recursive=true)
    end
    return destination
end

function _load_checkpoint(path, lattice; gauge, hz, Q, sites, expected_theta,
                          expected_settings, status)
    execution_identity = _execution_identity()
    _checkpoint_require(status in (:trial, :accepted), "unknown checkpoint status")
    _checkpoint_require(startswith(basename(normpath(path)), "checkpoint-") &&
        basename(dirname(normpath(path))) == string(status),
        "checkpoint status or published directory name mismatch")
    checksums = TOML.parsefile(joinpath(path, "checksums.toml"))
    _checkpoint_require(checksums["schema_version"] == 1, "unsupported checksum schema")
    for name in ("metadata.toml", "state.jls")
        file = joinpath(path, name)
        _checkpoint_require(filesize(file) == checksums["files"][name]["bytes"] &&
            _file_sha256(file) == checksums["files"][name]["sha256"],
            "checkpoint checksum or byte length mismatch: $name")
    end
    metadata = TOML.parsefile(joinpath(path, "metadata.toml"))
    _checkpoint_require(metadata["schema_version"] == 1 && metadata["format"] == _CHECKPOINT_FORMAT,
                        "unsupported checkpoint schema or format")
    _checkpoint_require(metadata["status"] == string(status), "checkpoint status mismatch")
    _checkpoint_require(metadata["runtime"] == Dict(execution_identity.runtime), "checkpoint runtime mismatch")
    _checkpoint_require(metadata["provenance"]["source_sha256"] ==
        Dict(execution_identity.source_sha256) &&
        metadata["provenance"]["manifest"] == execution_identity.manifest &&
        get(metadata["provenance"], "environment_sha256", nothing) ==
            Dict(execution_identity.environment_sha256),
        "checkpoint source or dependency manifest mismatch")
    charge = _checkpoint_charge(nsites(lattice), metadata["configuration"]["Q"]; Q)
    config = _checkpoint_configuration(lattice, gauge, hz; Q=charge)
    _checkpoint_require(metadata["configuration"] == config, "checkpoint model or geometry mismatch")
    N, Q = config["N"], config["Q"]
    state = metadata["state"]
    _checkpoint_require(_checkpoint_charge(N, state["Q"]) == Q,
                        "checkpoint charge metadata mismatch")
    theta = _finite_theta(state["theta"])
    expected_theta === nothing || _checkpoint_require(theta == expected_theta,
                                                     "checkpoint flux mismatch")
    theta_path = _checkpoint_flux_path(metadata["theta_path"], theta)
    settings = _checkpoint_settings(metadata["settings"])
    expected_settings === nothing || _checkpoint_require(
        settings == _checkpoint_settings(expected_settings), "checkpoint solver settings mismatch")
    diagnostics = _checkpoint_diagnostics(state, settings, N, Q)
    _checkpoint_require(metadata["baseline"]["theta"] == 0, "baseline must have zero flux")
    baseline_sz = _checkpoint_profile(metadata["baseline"]["sz"], N, Q)
    _checkpoint_require(length(theta_path) != 1 ||
        isapprox(diagnostics.sz, baseline_sz; atol=1e-10, rtol=0),
        "initial checkpoint baseline does not match its state")
    # Only trusted, integrity-checked, version-matched bytes reach deserialize.
    payload = open(deserialize, joinpath(path, "state.jls"))
    _checkpoint_require(length(payload.psi) == N, "checkpoint MPS length mismatch")
    _checkpoint_state_check(payload.psi, payload.sites, Q)
    _checkpoint_require(_checkpoint_site_record.(payload.sites) == state["site_indices"],
                        "checkpoint site identity metadata mismatch")
    sites === nothing || _checkpoint_require(_checkpoint_site_record.(sites) == state["site_indices"],
                                             "requested site indices do not match checkpoint")
    _checkpoint_require(isapprox(sz_profile(payload.psi), diagnostics.sz; atol=1e-10, rtol=0),
                        "checkpoint stored profile differs from MPS")
    _checkpoint_energy_check(payload.psi, payload.sites, lattice, theta,
                             gauge, config["hz"], diagnostics.energy)
    _checkpoint_require(length(payload.baseline_psi) == N, "baseline MPS length mismatch")
    _checkpoint_state_check(payload.baseline_psi, payload.sites, Q)
    _checkpoint_require(isapprox(sz_profile(payload.baseline_psi), baseline_sz; atol=1e-10, rtol=0),
                        "baseline profile does not match the zero-flux MPS")
    baseline = (; psi=payload.baseline_psi, theta=0.0, sz=baseline_sz,
                 sites=payload.sites, lattice=deepcopy(lattice),
                 hz=copy(config["hz"]), gauge, Q, execution_identity)
    return (; psi=payload.psi, sites=payload.sites, theta, gauge, Q,
             lattice=deepcopy(lattice), hz=copy(config["hz"]), settings,
             baseline, theta_path, diagnostics, metadata, execution_identity)
end

"""
    load_checkpoint(path, lattice; gauge=:seam, hz=nothing, Q=nothing, sites=nothing,
                    expected_theta=nothing, expected_settings=nothing, status=:accepted)

Validate and load a trusted local snapshot. The requested lattice (including
all oriented bonds), field and gauge must match exactly. Optional expectations
also compare unwrapped flux, solver settings and site identities. Check schema,
runtime, source/manifest hashes and file checksums before deserialization;
validate state normalization, charge, site indices and profile afterward.
Rebuild the stored Hamiltonian to verify energy, and check the saved zero-flux
MPS against its baseline profile. Omitted `Q` inherits the validated stored
integer charge; an explicit `Q` must be valid for the lattice and match it.
`status=:trial` permits explicit inspection of trials. Acceptance itself is
caller-supplied and does not establish convergence or branch continuity.
"""
function load_checkpoint(path::AbstractString, lattice::KagomeCylinder;
                         gauge::Symbol=:seam, hz=nothing, Q=nothing, sites=nothing,
                         expected_theta=nothing, expected_settings=nothing,
                         status::Symbol=:accepted)
    try
        return _load_checkpoint(path, lattice; gauge, hz, Q, sites, expected_theta,
                                expected_settings, status)
    catch exception
        exception isa InterruptException && rethrow()
        exception isa ArgumentError && rethrow()
        throw(ArgumentError("invalid checkpoint ($(nameof(typeof(exception))))"))
    end
end

"""
    resume_dmrg(path, lattice, theta; gauge=:seam, hz=nothing, Q=nothing,
                expected_settings=nothing, outputlevel=0)

Start a fresh DMRG batch from an accepted checkpoint using its saved solver
settings, validated charge and exact site indices. An explicit `Q` must match
the checkpoint; omission inherits it. Rebuild the Hamiltonian and environments at
the requested unwrapped `theta`. This is a restart helper, not a branch tracker
or a continuation acceptance rule. It does not overwrite the checkpoint.
"""
function resume_dmrg(path::AbstractString, lattice::KagomeCylinder, theta::Real;
                     gauge::Symbol=:seam, hz=nothing, Q=nothing, expected_settings=nothing,
                     outputlevel::Integer=0)
    saved = load_checkpoint(path, lattice; gauge, hz, Q, expected_settings)
    s = saved.settings
    return run_dmrg(lattice, theta; sites=saved.sites, psi0=saved.psi, Q=saved.Q,
        seed=s.seed, initial_linkdim=s.initial_linkdim, nsweeps=s.nsweeps,
        maxdim=s.maxdim, cutoff=s.cutoff, noise=s.noise, eigsolve_tol=s.eigsolve_tol,
        eigsolve_krylovdim=s.eigsolve_krylovdim, eigsolve_maxiter=s.eigsolve_maxiter,
        measure_variance=s.measure_variance, gauge, hz, outputlevel)
end
