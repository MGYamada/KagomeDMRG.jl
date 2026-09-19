#!/usr/bin/env julia
# Bounded NN-sector API study, separate from the package regression suite.
# julia --project=. --startup-file=no --threads=1 examples/validate_charge_sectors.jl NEW_OUTPUT
using KagomeDMRG, ITensors, ITensorMPS, LinearAlgebra, SHA, TOML, Dates
BLAS.set_num_threads(1)
ITensors.disable_threaded_blocksparse()
const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "test", "reference_ed.jl"))
include(joinpath(ROOT, "test", "itensor_helpers.jl"))

length(ARGS) == 1 || error("usage: validate_charge_sectors.jl new-output-directory")
const OUTPUT = abspath(only(ARGS))
ispath(OUTPUT) && error("choose a new output directory")
mkpath(OUTPUT)
const CASES = ((Lx=1, Ly=4, Q=0, theta=0.0),
               (Lx=1, Ly=4, Q=0, theta=0.37),
               (Lx=1, Ly=3, Q=-1, theta=0.0),
               (Lx=1, Ly=3, Q=1, theta=0.0),
               (Lx=1, Ly=3, Q=3, theta=0.0))
const SOLVER = (; seed=11, nsweeps=8, maxdim=64, cutoff=0.0, noise=0.0,
    eigsolve_tol=1e-12, eigsolve_krylovdim=30, eigsolve_maxiter=10,
    measure_variance=true)
const LIMITS = Dict("energy_error_per_site"=>1e-8, "residual_norm"=>2e-6,
    "ground_space_leakage"=>2e-6, "max_sz_error"=>1e-6,
    "total_sz_error"=>1e-12, "state_norm_error"=>1e-12,
    "gauge_matrix_error"=>1e-11, "mpo_action_error"=>1e-11,
    "uniform_mpo_action_error"=>1e-11)
const EXTRA_SOURCES = ("examples/validate_charge_sectors.jl", "test/reference_ed.jl",
                       "test/itensor_helpers.jl")
extra_hashes() = Dict(p=>bytes2hex(sha256(read(joinpath(ROOT, p)))) for p in EXTRA_SOURCES)
const EXTRA_BEFORE = extra_hashes()
const START = time_ns()
elapsed() = (time_ns()-START) / 1e9
const RECORD = Dict{String,Any}("schema_version"=>1, "status"=>"running",
    "recorded_at_utc"=>string(now(UTC)), "scope"=>"explicit_charge_small_system_validation",
    "model"=>"nearest_neighbor_isotropic_Heisenberg_J1=1_hz=0",
    "phase_identification"=>"not_attempted", "plateau_identification"=>"not_attempted",
    "axial_pump"=>"not_measured_Lx=1", "trajectory_acceptance"=>"not_attempted",
    "code"=>KagomeDMRG._checkpoint_provenance(), "extra_source_sha256"=>EXTRA_BEFORE,
    "runtime"=>merge(KagomeDMRG._checkpoint_runtime(), Dict(
        "julia_threads"=>Threads.nthreads(), "blas_threads"=>BLAS.get_num_threads())),
    "solver"=>Dict(string(k)=>v for (k,v) in pairs(SOLVER)), "limits"=>LIMITS,
    "planned_cases"=>[Dict(string(k)=>v for (k,v) in pairs(c)) for c in CASES],
    "time_limit_seconds"=>300, "time_limit_enforcement"=>"between_bounded_points",
    "runs"=>Dict{String,Any}[])
persist() = KagomeDMRG._write_flux_record(joinpath(OUTPUT, "validation.toml"), RECORD)

# Contract H*probe directly, retaining only one physical index per site.
# This avoids materializing the full 4096-by-4096 twelve-site MPO matrix.
function mpo_sector_action(H, probe, sites, basis)
    length(sites) <= 12 || error("this study bounds MPO action to twelve sites")
    tensor = _contract_small_system([H[i]*probe[i] for i in eachindex(sites)])
    amplitudes = vec(Array(tensor, prime.(sites)...))
    return amplitudes[2^length(sites) .- Int.(basis)]
end

function measure_case(case, baseline)
    (; Lx, Ly, Q, theta) = case
    N = 3Lx*Ly
    lattice = kagome_cylinder(Lx, Ly)
    ed_seconds = @elapsed begin
        ref = reference_hamiltonian(Lx, Ly, theta; Q)
        eig = eigen(Hermitian(ref.H))
    end
    dmrg_seconds = @elapsed result = baseline === nothing ?
        run_dmrg(lattice, theta; Q, SOLVER...) :
        run_dmrg(lattice, theta; Q, sites=baseline.sites, psi0=baseline.psi, SOLVER...)
    cluster = findall(e -> e-first(eig.values) < 1e-10, eig.values)
    ground = eig.vectors[:, cluster]
    v = sector_amplitudes(result.psi, result.sites, ref.basis)
    projection = ground*(ground'*v)
    projected_norm = norm(projection)
    projected_norm > 0 || error("zero ground-space projection")
    matched = projection/projected_norm
    uniform = reference_hamiltonian(Lx, Ly, theta; Q, gauge=:uniform)
    rotation = Diagonal(reference_gauge_diagonal(ref.basis, Lx, Ly, theta))
    probe = initial_mps(result.sites; Q, seed=33)
    probe_v = sector_amplitudes(probe, result.sites, ref.basis)
    uniform_H = twisted_exchange_mpo(result.sites, lattice, theta; gauge=:uniform)
    row = Dict{String,Any}(
        "configuration"=>KagomeDMRG._checkpoint_configuration(lattice, :seam, nothing; Q),
        "theta"=>theta, "Q"=>Q, "N"=>N, "basis_dimension"=>length(ref.basis),
        "dmrg_seconds"=>dmrg_seconds, "ed_seconds"=>ed_seconds,
        "ed_energy"=>first(eig.values), "energy"=>result.energy,
        "ground_degeneracy"=>length(cluster),
        "gap_above_ground_space"=>eig.values[length(cluster)+1]-first(eig.values),
        "ed_ground_residual_norm"=>norm(ref.H*ground-ground*Diagonal(eig.values[cluster])),
        "energy_error_per_site"=>abs(result.energy-first(eig.values))/N,
        "residual_norm"=>norm(ref.H*v-result.energy*v),
        "ground_space_leakage"=>norm(v-projection),
        "max_sz_error"=>maximum(abs.(result.sz-reference_sz(matched, ref.basis, N))),
        "total_sz_error"=>abs(sum(result.sz)-Q/2), "state_norm_error"=>abs(norm(v)-1),
        "gauge_matrix_error"=>norm(uniform.H-rotation*ref.H*rotation'),
        "mpo_action_error"=>norm(mpo_sector_action(result.H, probe, result.sites, ref.basis)-ref.H*probe_v),
        "uniform_mpo_action_error"=>norm(mpo_sector_action(uniform_H, probe, result.sites, ref.basis)-uniform.H*probe_v),
        "sz_profile"=>result.sz, "raw_variance"=>result.variance,
        "maxlinkdim"=>maxlinkdim(result.psi), "sweep_energies"=>result.sweep_energies,
        "measured_truncation_errors"=>result.max_truncation_errors)
    row["passed"] = all(isfinite(row[k]) && row[k] <= limit for (k,limit) in LIMITS) &&
                    row["ed_ground_residual_norm"] < 1e-10
    zero = baseline === nothing ? result : baseline
    path = save_checkpoint(OUTPUT, result; baseline=zero,
        theta_path=theta == 0 ? [0.0] : [0.0, theta], status=:trial)
    loaded = load_checkpoint(path, lattice; status=:trial, expected_theta=theta)
    row["checkpoint"] = relpath(path, OUTPUT)
    row["checkpoint_inherited_Q"] = loaded.Q
    row["checkpoint_overlap_error"] = abs(abs(inner(loaded.psi, result.psi))-1)
    row["passed"] &= loaded.Q == Q && row["checkpoint_overlap_error"] < 1e-12
    return result, row
end

function main()
    baseline = nothing
    persist()
    for case in CASES
        if elapsed() > RECORD["time_limit_seconds"]
            RECORD["status"] = "time_budget_exceeded"
            return
        end
        try
            result, row = measure_case(case, case.theta == 0 ? nothing : baseline)
            case.Ly == 4 && case.theta == 0 && (baseline = result)
            push!(RECORD["runs"], row)
            persist()
            println("N=", row["N"], " Q=", row["Q"], " theta=", row["theta"],
                " energy/site error=", row["energy_error_per_site"],
                " residual=", row["residual_norm"], " passed=", row["passed"])
        catch exception
            exception isa InterruptException && rethrow()
            exception isa OutOfMemoryError && rethrow()
            push!(RECORD["runs"], Dict("case"=>Dict(string(k)=>v for (k,v) in pairs(case)),
                                      "status"=>"failed", "error_type"=>string(nameof(typeof(exception)))))
            RECORD["status"] = "failed"
            rethrow()
        end
    end
    neighbor = Dict(r["Q"]=>r for r in RECORD["runs"] if r["N"] == 9)
    RECORD["neighboring_sector_comparison"] = Dict(
        "scope"=>"finite_N9_adjacent_sector_energy_differences_only",
        "spin_reversal_energy_error"=>abs(neighbor[1]["ed_energy"]-neighbor[-1]["ed_energy"]),
        "h_minus"=>neighbor[1]["ed_energy"]-neighbor[-1]["ed_energy"],
        "h_plus"=>neighbor[3]["ed_energy"]-neighbor[1]["ed_energy"])
    RECORD["status"] = all(r["passed"] for r in RECORD["runs"]) ? "passed" : "outside_limits"
end

try
    main()
catch exception
    RECORD["status"] = exception isa InterruptException ? "interrupted" : "failed"
    RECORD["exception_type"] = string(nameof(typeof(exception)))
    rethrow()
finally
    RECORD["elapsed_seconds"] = elapsed()
    RECORD["extra_sources_unchanged"] = extra_hashes() == EXTRA_BEFORE
    RECORD["extra_sources_unchanged"] || (RECORD["status"] = "source_changed")
    persist()
end
RECORD["status"] == "passed" || error("charge-sector study did not pass; inspect validation.toml")
