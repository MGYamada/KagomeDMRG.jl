#!/usr/bin/env julia
# Two bounded single-point optimizations, not a CSL phase or pump experiment.
# julia --project=. --startup-file=no --threads=1 examples/validate_extended_model.jl NEW_OUTPUT
using KagomeDMRG, ITensors, ITensorMPS, LinearAlgebra, SHA, TOML, Dates
BLAS.set_num_threads(1)
ITensors.disable_threaded_blocksparse()
const ROOT = normpath(joinpath(@__DIR__, ".."))
include(joinpath(ROOT, "test", "reference_ed.jl"))
include(joinpath(ROOT, "test", "reference_extended_ed.jl"))
include(joinpath(ROOT, "test", "itensor_helpers.jl"))

length(ARGS) == 1 || error("usage: validate_extended_model.jl new-output-directory")
const OUTPUT = abspath(only(ARGS))
ispath(OUTPUT) && error("choose a new output directory")
mkpath(OUTPUT)
const COUPLINGS = (; J1=1.0, J2=0.5, J3=0.5)
const SOLVER = (; seed=11, nsweeps=8, maxdim=64, cutoff=0.0, noise=0.0,
    eigsolve_tol=1e-12, eigsolve_krylovdim=30, eigsolve_maxiter=10,
    measure_variance=true)
const LIMITS = Dict("energy_error_per_site"=>1e-8, "residual_norm"=>2e-6,
    "ground_space_leakage"=>2e-6, "max_sz_error"=>1e-6, "max_bond_energy_error"=>1e-6,
    "total_sz_error"=>1e-12, "state_norm_error"=>1e-12,
    "gauge_matrix_error"=>1e-11, "mpo_action_error"=>1e-11,
    "uniform_mpo_action_error"=>1e-11, "uniform_state_error"=>1e-11,
    "bond_energy_gauge_error"=>1e-11, "bond_energy_sum_error"=>1e-11,
    "nn_reduction_action_error"=>1e-11, "periodicity_action_error"=>1e-11)
const EXTRA_SOURCES = ("examples/validate_extended_model.jl", "test/reference_ed.jl",
    "test/reference_extended_ed.jl", "test/itensor_helpers.jl")
extra_hashes() = Dict(p=>bytes2hex(sha256(read(joinpath(ROOT, p)))) for p in EXTRA_SOURCES)
const EXTRA_BEFORE = extra_hashes()
const START = time_ns()
elapsed() = (time_ns()-START)/1e9
const RECORD = Dict{String,Any}("schema_version"=>1, "status"=>"running",
    "recorded_at_utc"=>string(now(UTC)), "scope"=>"extended_exchange_small_system_validation",
    "model"=>"J1_J2_J3_hexagon_opposite_only_no_field_or_chirality_seed",
    "couplings"=>Dict(string(k)=>v for (k,v) in pairs(COUPLINGS)),
    "Lx"=>1, "Ly"=>4, "N"=>12, "Q"=>0, "planned_theta"=>[0.0, 0.37],
    "phase_identification"=>"not_attempted", "axial_pump"=>"not_measured_Lx=1",
    "trajectory_acceptance"=>"not_attempted", "chirality"=>"not_measured",
    "code"=>KagomeDMRG._checkpoint_provenance(), "extra_source_sha256"=>EXTRA_BEFORE,
    "runtime"=>merge(KagomeDMRG._checkpoint_runtime(), Dict(
        "julia_threads"=>Threads.nthreads(), "blas_threads"=>BLAS.get_num_threads())),
    "solver"=>Dict(string(k)=>v for (k,v) in pairs(SOLVER)), "limits"=>LIMITS,
    "time_limit_seconds"=>300, "time_limit_enforcement"=>"between_bounded_points",
    "runs"=>Dict{String,Any}[])
persist() = KagomeDMRG._write_flux_record(joinpath(OUTPUT, "validation.toml"), RECORD)

function mpo_action(H, psi, sites, basis)
    length(sites) == 12 || error("this study is restricted to twelve sites")
    tensor = _contract_small_system([H[i]*psi[i] for i in eachindex(sites)])
    return vec(Array(tensor, prime.(sites)...))[2^length(sites) .- Int.(basis)]
end

bond_key(b) = b.i < b.j ? (b.i, b.j, b.wy) : (b.j, b.i, -b.wy)
function reference_bond_energies(v, basis, bonds, theta)
    return Dict(bond_key(b) => b.Jz*real(reference_correlation(v,basis,12,b.i,b.j)) +
        b.Jxy*real(cis(b.wy*theta)*reference_correlation(v,basis,12,b.i,b.j;
                                                       operators=("S+", "S-"))) for b in bonds)
end

function measure_point(theta, baseline)
    lattice = kagome_j1j2j3_cylinder(1, 4; COUPLINGS...)
    ed_seconds = @elapsed begin
        ref = reference_extended_hamiltonian(1, 4, theta; Q=0, COUPLINGS...)
        eig = eigen(Hermitian(ref.H))
    end
    dmrg_seconds = @elapsed result = baseline === nothing ?
        run_dmrg(lattice, theta; Q=0, SOLVER...) :
        run_dmrg(lattice, theta; Q=0, sites=baseline.sites, psi0=baseline.psi, SOLVER...)
    v = sector_amplitudes(result.psi, result.sites, ref.basis)
    cluster = findall(e -> e-first(eig.values) < 1e-10, eig.values)
    ground = eig.vectors[:,cluster]
    projection = ground*(ground'*v)
    norm(projection) > 0 || error("zero ground-space projection")
    matched = projection/norm(projection)
    independent_bonds = reference_extended_bonds(1, 4; COUPLINGS...)
    expected_bonds = reference_bond_energies(matched, ref.basis, independent_bonds, theta)
    energies = bond_energies(result.psi, lattice, theta)
    families = bond_families(lattice)
    uniform = reference_extended_hamiltonian(1, 4, theta; Q=0, gauge=:uniform, COUPLINGS...)
    rotation = Diagonal(reference_gauge_diagonal(ref.basis, 1, 4, theta))
    uniform_psi = deepcopy(result.psi)
    for (i, chi) in enumerate(gauge_angles(lattice, theta))
        s = result.sites[i]
        U = cos(chi/2)*op("Id",s) + 2im*sin(chi/2)*op("Sz",s)
        uniform_psi[i] = noprime(U*uniform_psi[i])
    end
    probe = initial_mps(result.sites; Q=0, seed=33)
    probe_v = sector_amplitudes(probe, result.sites, ref.basis)
    uniform_H = twisted_exchange_mpo(result.sites, lattice, theta; gauge=:uniform)
    periodic_H = twisted_exchange_mpo(result.sites, lattice, theta+2pi)
    nn = kagome_cylinder(1, 4)
    zero_extended = kagome_j1j2j3_cylinder(1, 4; J2=0, J3=0)
    nn_action = mpo_action(twisted_exchange_mpo(result.sites, nn, theta), probe, result.sites, ref.basis)
    zero_action = mpo_action(twisted_exchange_mpo(result.sites, zero_extended, theta), probe, result.sites, ref.basis)
    seam_action = mpo_action(result.H, probe, result.sites, ref.basis)
    row = Dict{String,Any}(
        "configuration"=>KagomeDMRG._checkpoint_configuration(lattice,:seam,nothing;Q=0),
        "theta"=>theta, "basis_dimension"=>length(ref.basis),
        "dmrg_seconds"=>dmrg_seconds, "ed_seconds"=>ed_seconds,
        "energy"=>result.energy, "ed_energy"=>first(eig.values),
        "ground_degeneracy"=>length(cluster),
        "gap_above_ground_space"=>eig.values[length(cluster)+1]-first(eig.values),
        "ed_ground_residual_norm"=>norm(ref.H*ground-ground*Diagonal(eig.values[cluster])),
        "energy_error_per_site"=>abs(result.energy-first(eig.values))/12,
        "residual_norm"=>norm(ref.H*v-result.energy*v),
        "ground_space_leakage"=>norm(v-projection),
        "max_sz_error"=>maximum(abs.(result.sz-reference_sz(matched,ref.basis,12))),
        "max_bond_energy_error"=>maximum(abs(energies[k]-expected_bonds[bond_key(b)])
                                         for (k,b) in enumerate(lattice.bonds)),
        "total_sz_error"=>abs(sum(result.sz)), "state_norm_error"=>abs(norm(v)-1),
        "gauge_matrix_error"=>norm(uniform.H-rotation*ref.H*rotation'),
        "mpo_action_error"=>norm(seam_action-ref.H*probe_v),
        "uniform_mpo_action_error"=>norm(mpo_action(uniform_H,probe,result.sites,ref.basis)-uniform.H*probe_v),
        "uniform_state_error"=>norm(sector_amplitudes(uniform_psi,result.sites,ref.basis)-rotation*v),
        "bond_energy_gauge_error"=>maximum(abs.(bond_energies(uniform_psi,lattice,theta;gauge=:uniform)-energies)),
        "bond_energy_sum_error"=>abs(sum(energies)-result.energy),
        "nn_reduction_action_error"=>norm(nn_action-zero_action),
        "periodicity_action_error"=>norm(mpo_action(periodic_H,probe,result.sites,ref.basis)-seam_action),
        "sz_profile"=>result.sz, "bond_energy_profile"=>energies,
        "family_energy"=>Dict(string(f)=>sum(energies[families .== f]) for f in (:J1,:J2,:J3)),
        "raw_variance"=>result.variance, "maxlinkdim"=>maxlinkdim(result.psi),
        "sweep_energies"=>result.sweep_energies,
        "measured_truncation_errors"=>result.max_truncation_errors)
    snapshot = save_checkpoint(OUTPUT,result; baseline=baseline === nothing ? result : baseline,
        theta_path=theta == 0 ? [0.0] : [0.0,theta],status=:trial)
    loaded = load_checkpoint(snapshot,lattice;status=:trial,expected_theta=theta)
    row["checkpoint"] = relpath(snapshot,OUTPUT)
    row["checkpoint_overlap_error"] = abs(abs(inner(loaded.psi,result.psi))-1)
    row["checkpoint_families_match"] = getindex.(loaded.metadata["configuration"]["bonds"],"family") == string.(families)
    row["nn_checkpoint_request_rejected"] = try
        load_checkpoint(snapshot,nn;Q=0,status=:trial)
        false
    catch exception
        exception isa ArgumentError || rethrow()
        true
    end
    row["passed"] = all(isfinite(row[k]) && row[k] <= limit for (k,limit) in LIMITS) &&
        row["ed_ground_residual_norm"] < 1e-10 && loaded.Q == 0 &&
        row["checkpoint_overlap_error"] < 1e-12 && row["checkpoint_families_match"] &&
        row["nn_checkpoint_request_rejected"]
    return result,row
end

function main()
    baseline = nothing
    persist()
    for theta in (0.0,0.37)
        if elapsed() > RECORD["time_limit_seconds"]
            RECORD["status"] = "time_budget_exceeded"
            return
        end
        result,row = measure_point(theta,baseline)
        baseline === nothing && (baseline = result)
        push!(RECORD["runs"],row)
        persist()
        println("theta=",theta," residual=",row["residual_norm"],
                " bond error=",row["max_bond_energy_error"]," passed=",row["passed"])
    end
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
RECORD["status"] == "passed" || error("extended-model study did not pass; inspect validation.toml")
