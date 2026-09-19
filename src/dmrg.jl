"""
    initial_mps(sites; Q=nothing, seed=0, linkdim=4)

Make a reproducible complex random MPS with integer total charge `Q=2Sz`.
Omitting `Q` selects `M/Msat=1/9`; see [`target_sector`](@ref).
A private random-number generator chooses a product configuration with
exactly `(N+Q)/2` up spins before randomizing within the same total charge.
"""
function initial_mps(sites; Q=nothing, seed::Integer=0, linkdim::Integer=4)
    sector = target_sector(length(sites); Q)
    _check_sites(sites, sector.N)
    linkdim > 0 || throw(ArgumentError("linkdim must be positive"))
    rng = MersenneTwister(seed)
    labels = fill("Dn", sector.N)
    labels[randperm(rng, sector.N)[1:sector.Nup]] .= "Up"
    psi = random_mps(rng, ComplexF64, sites, labels; linkdims=linkdim)
    flux(psi) == QN("Sz", sector.Q) || error("initial MPS has incorrect total charge")
    return psi
end

# Collect measured factorization truncation errors on both halves of each sweep.
# The configured cutoff is kept separately in the returned settings.
mutable struct _SweepDiagnostics{F} <: AbstractObserver
    energies::Vector{Float64}
    max_truncation_errors::Vector{Float64}
    progress_callback::F
    started_ns::UInt64
end
_SweepDiagnostics(callback=nothing, started_ns::UInt64=time_ns()) =
    _SweepDiagnostics(Float64[], Float64[], callback, started_ns)

function _dmrg_progress(callback, started_ns; event...)
    callback === nothing && return nothing
    callback((; event..., elapsed_seconds=Float64(time_ns() - started_ns) * 1e-9))
    return nothing
end

function ITensorMPS.measure!(obs::_SweepDiagnostics; sweep, half_sweep, bond,
                            energy, spec, psi, kwargs...)
    # The orthogonality center follows the sweep. Some upstream QN
    # factorizations can remove every block at a degenerate cutoff, while
    # still reporting a finite truncation error. Never accept that state.
    center = half_sweep == 1 ? bond + 1 : bond
    center_norm = norm(psi[center])
    isfinite(center_norm) && center_norm > 0 || error(
        "DMRG factorization produced a zero or nonfinite state at sweep $sweep, " *
        "half-sweep $half_sweep, bond $bond. Increase maxdim and check the " *
        "upstream QN truncation boundary before continuing.")
    # Equal or nearly equal weights at a cross-sector cutoff may cause the
    # backend to retain fewer states than its Spectrum reports. A nonzero
    # norm alone does not detect that underreported discarded probability.
    spectrum = eigs(spec)
    spectrum !== nothing && length(spectrum) == dim(linkind(psi, bond)) || error(
        "DMRG truncation spectrum does not match the retained bond dimension " *
        "at sweep $sweep, half-sweep $half_sweep, bond $bond. Increase maxdim " *
        "to retain the degenerate boundary before continuing.")
    while length(obs.max_truncation_errors) < sweep
        push!(obs.max_truncation_errors, 0.0)
    end
    obs.max_truncation_errors[sweep] =
        max(obs.max_truncation_errors[sweep], truncerror(spec))
    obs.progress_callback === nothing || _dmrg_progress(
        obs.progress_callback, obs.started_ns; kind=:bond,
        sweep=Int(sweep), half_sweep=Int(half_sweep), bond=Int(bond),
        energy=Float64(real(energy)), truncation_error=Float64(truncerror(spec)),
        retained_dimension=dim(linkind(psi, bond)))
    if half_sweep == 2 && bond == 1
        push!(obs.energies, real(energy))
        obs.progress_callback === nothing || _dmrg_progress(
            obs.progress_callback, obs.started_ns; kind=:sweep, sweep=Int(sweep),
            energy=Float64(real(energy)),
            max_truncation_error=obs.max_truncation_errors[sweep], maxlinkdim=maxlinkdim(psi))
    end
    return nothing
end

"""
    run_dmrg(lattice, theta; Q=nothing, sites=spin_sites(lattice), psi0=nothing, ...)

Run the ITensor two-site U(1) reference solver at a single, unwrapped flux.
The result contains the optimized `psi`, freshly built `H`, energy, `sz`,
per-sweep measured truncation errors, raw energy variance, and solver settings.
`psi0` is copied; its site indices and total charge must match the request.
Omitting `Q` selects the 1/9 sector, including when `psi0` is supplied.
Set `Q` explicitly for other sectors; the solver never infers it from `psi0`.

This is a single-point optimizer, not an adiabatic branch-tracking driver.
No convergence or phase label is inferred from a low energy or a warm start.
The raw variance may be slightly negative from floating-point cancellation.
Set `measure_variance=false` to skip the additional H² contraction explicitly.
Local Krylov convergence flags are not exposed by the upstream `dmrg` API.
With nonzero `noise`, truncation errors refer to a perturbed density matrix;
they must not be interpreted as exact discarded wavefunction probabilities.
An empty or nonfinite state after factorization raises an error immediately.
The pinned local NDTensors correction allocates retained ranks from one global
selection, including ties across QN sectors. A hard `maxdim` can split a tied
subspace; equal weights use block-coordinate then local-index order. Guards
still reject zero states and reported-spectrum/retained-dimension mismatches.
The calibrated truncation loss is not an error bound on physical observables.

`progress_callback(event)` optionally receives scalar-only NamedTuples with
`kind=:phase` (initialization, mpo, sweeps, energy, variance, sz; started/completed),
`:bond` (after each guarded local update), or `:sweep` (after each full sweep).
`elapsed_seconds` is cumulative wall time since entry to this function, including
callback time. Bond/sweep energies are local Ritz values, not final expectations.
Skipped variance has no phase events. Callback exceptions propagate; callbacks
receive no solver state and are absent from returned settings and checkpoints.
"""
function run_dmrg(lattice::KagomeCylinder, theta::Real;
                  Q=nothing, sites=spin_sites(lattice), psi0=nothing, seed::Integer=0,
                  initial_linkdim::Integer=4, nsweeps::Integer=8,
                  maxdim=[16, 32, 64], cutoff::Real=1e-12, noise::Real=0.0,
                  eigsolve_tol::Real=1e-12, eigsolve_krylovdim::Integer=20,
                  eigsolve_maxiter::Integer=10, gauge::Symbol=:seam, hz=nothing,
                  outputlevel::Integer=0, measure_variance::Bool=true,
                  progress_callback=nothing)
    started_ns = time_ns()
    phase(name, status) = _dmrg_progress(progress_callback, started_ns;
        kind=:phase, phase=name, status)
    phase(:initialization, :started)
    execution_identity = _execution_identity()
    sector = target_sector(nsites(lattice); Q)
    _check_sites(sites, sector.N)
    nsweeps > 0 || throw(ArgumentError("nsweeps must be positive"))
    dims = maxdim isa Integer ? [Int(maxdim)] : Int.(collect(maxdim))
    !isempty(dims) && all(>(0), dims) || throw(ArgumentError("maxdim must be positive"))
    isfinite(cutoff) && cutoff >= 0 || throw(ArgumentError("cutoff must be finite and nonnegative"))
    isfinite(noise) && noise >= 0 || throw(ArgumentError("noise must be finite and nonnegative"))
    isfinite(eigsolve_tol) && eigsolve_tol > 0 || throw(ArgumentError("eigsolve_tol must be positive"))
    eigsolve_krylovdim >= 2 || throw(ArgumentError("eigsolve_krylovdim must be at least 2"))
    eigsolve_maxiter > 0 || throw(ArgumentError("eigsolve_maxiter must be positive"))
    if psi0 === nothing
        psi = initial_mps(sites; Q=sector.Q, seed, linkdim=initial_linkdim)
    else
        length(psi0) == sector.N || throw(ArgumentError("initial MPS has incorrect length"))
        all(i -> siteind(psi0, i) == sites[i], eachindex(sites)) ||
            throw(ArgumentError("initial MPS site indices do not match"))
        flux(psi0) == QN("Sz", sector.Q) ||
            throw(ArgumentError("initial MPS must have the requested integer charge Q=$(sector.Q)"))
        psi = complex(deepcopy(psi0))
    end
    # Keep the model attached to the result for checkpoint validation. A
    # caller mutating its lattice or field vector later must not relabel it.
    saved_lattice = deepcopy(lattice)
    fields = hz === nothing ? zeros(sector.N) : Float64.(hz)
    phase(:initialization, :completed)
    phase(:mpo, :started)
    H = twisted_exchange_mpo(sites, saved_lattice, theta; gauge, hz=fields)
    phase(:mpo, :completed)
    obs = _SweepDiagnostics(progress_callback, started_ns)
    phase(:sweeps, :started)
    local_energy, psi = dmrg(H, psi; nsweeps, maxdim=dims, cutoff, noise,
                            eigsolve_tol, eigsolve_krylovdim, eigsolve_maxiter,
                            observer=obs, outputlevel, ishermitian=true)
    normalize!(psi)
    flux(psi) == QN("Sz", sector.Q) || error("DMRG changed the total charge")
    phase(:sweeps, :completed)
    phase(:energy, :started)
    energy_complex = inner(psi', H, psi)
    abs(imag(energy_complex)) < 1e-10 * max(1, abs(real(energy_complex))) ||
        error("energy has an unexpectedly large imaginary part")
    energy = real(energy_complex)
    phase(:energy, :completed)
    variance = nothing
    if measure_variance
        phase(:variance, :started)
        variance = real(inner(H, psi, H, psi)) - energy^2
        phase(:variance, :completed)
    end
    phase(:sz, :started)
    sz = sz_profile(psi)
    phase(:sz, :completed)
    settings = (; seed=Int(seed), initial_linkdim=Int(initial_linkdim),
                 nsweeps=Int(nsweeps), maxdim=dims, cutoff=Float64(cutoff),
                 noise=Float64(noise), eigsolve_tol=Float64(eigsolve_tol),
                 eigsolve_krylovdim=Int(eigsolve_krylovdim),
                 eigsolve_maxiter=Int(eigsolve_maxiter), measure_variance,
                 initialization=psi0 === nothing ? "random_fixed_charge" : "provided_mps")
    execution_identity == _execution_identity() ||
        error("implementation changed during DMRG")
    return (; psi, H, sites, lattice=saved_lattice, hz=copy(fields),
             theta=Float64(theta), gauge, Q=sector.Q, energy,
             local_energy=real(local_energy), variance, sz,
             sweep_energies=obs.energies,
             max_truncation_errors=obs.max_truncation_errors, settings,
             execution_identity)
end
