#!/usr/bin/env julia
# Focused preparation/removal validation, not a VBC phase-identification test.
# julia --project=research --startup-file=no --threads=1 examples/validate_vbc_preparation.jl NEW_OUTPUT
const VBC_CHECK_START = time_ns()
using KagomeDMRG, ITensors, ITensorMPS, LinearAlgebra, SHA, TOML, Dates
BLAS.set_num_threads(1)
ITensors.disable_threaded_blocksparse()
const VBC_CHECK_ROOT = normpath(joinpath(@__DIR__, ".."))
include("vbc_motif_tools.jl")
include(joinpath(VBC_CHECK_ROOT, "test/reference_ed.jl"))
include(joinpath(VBC_CHECK_ROOT, "test/itensor_helpers.jl"))

vbc_check_key(b) = b.i < b.j ? (b.i, b.j, b.wy) : (b.j, b.i, -b.wy)
vbc_check_hash(path) = bytes2hex(sha256(read(path)))

"""Independent Cartesian-neighbor and bit-basis exchange construction.

Only the prescribed motif weights are shared input. This checks applying the
preparation to the Hamiltonian, its orientation/phases, and its exact removal;
the mapping from the literature drawing to these weights is a separate audit.
"""
function vbc_reference_preparation(lattice, weights, lambda, theta)
    N = nsites(lattice)
    N == 9 || error("this dense preparation reference is bounded to N9")
    bonds = reference_bonds(lattice.Lx, lattice.Ly)
    prescribed = Dict(vbc_check_key(b)=>w for (b,w) in zip(lattice.bonds, weights))
    Set(keys(prescribed)) == Set(vbc_check_key.(bonds)) || error("NN geometry differs from Cartesian reference")
    basis = _reference_basis(N, :all)
    H = zeros(ComplexF64, length(basis), length(basis))
    for (column, state) in enumerate(basis), b in bonds
        J = 1 + lambda * prescribed[vbc_check_key(b)]
        zi, zj = _reference_sz_bit(state,b.i), _reference_sz_bit(state,b.j)
        H[column,column] += J * zi * zj
        if zi != zj
            flipped = state ⊻ (UInt64(1) << (b.i-1)) ⊻ (UInt64(1) << (b.j-1))
            H[Int(flipped)+1,column] += (J/2) * cis((zi < zj ? 1 : -1) * b.wy * theta)
        end
    end
    return (; H, basis)
end

function vbc_preparation_checks(output)
    ispath(output) && !(readdir(output) in (String[],["execution.toml"])) &&
        error("choose a new, empty, or freshly supervised output directory")
    mkpath(output)
    isfile(joinpath(output,"execution.toml")) &&
        TOML.parsefile(joinpath(output,"execution.toml"))["status"]!="running" &&
        error("supervised validation needs a running execution record")
    Threads.nthreads()==BLAS.get_num_threads()==1 || error("one Julia/BLAS thread required")
    sources = ("examples/validate_vbc_preparation.jl", "examples/vbc_motif_tools.jl",
        "test/reference_ed.jl", "test/itensor_helpers.jl")
    hashes = Dict(p=>vbc_check_hash(joinpath(VBC_CHECK_ROOT,p)) for p in sources)
    record = Dict{String,Any}("schema_version"=>1,"status"=>"running",
        "recorded_at_utc"=>string(now(UTC)),"N"=>9,"Q"=>1,
        "matrix_theta"=>0.371,"warm_start_theta"=>0.0,"lambda"=>0.3,
        "analysis_source_sha256"=>hashes,"code"=>KagomeDMRG._checkpoint_provenance(),
        "runtime"=>KagomeDMRG._checkpoint_runtime(),
        "scope"=>"new_bond_preparation_and_removal_only_not_published_wavefunction_or_phase",
        "reference_independence"=>"Cartesian_neighbor_search_and_spin_bits;motif_weights_are_shared_input",
        "checks"=>Dict{String,Any}(),"metrics"=>Dict{String,Any}())
    persist() = KagomeDMRG._write_flux_record(joinpath(output,"validation.toml"),record)
    check(name,value) = (record["checks"][name]=Bool(value); value || error("failed preparation check: $name"))
    persist()
    for p in sources
        destination = joinpath(output,"analysis-sources",p)
        mkpath(dirname(destination))
        cp(joinpath(VBC_CHECK_ROOT,p),destination)
    end
    try
        lattice = kagome_cylinder(1,3)
        sites = spin_sites(lattice)
        theta, lambda = 0.371, 0.3
        reference_nn = reference_hamiltonian(1,3,theta;nup=:all)
        indices = 2^9 .- Int.(reference_nn.basis)
        prep_for_solve = nothing
        for kind in (:hourglass,:windmill)
            name = string(kind)
            template = vbc_bond_template(lattice;kind,origin=(0,0))
            prepared = vbc_prepared_lattice(lattice,template,lambda)
            removed = vbc_prepared_lattice(lattice,template,0.0)
            reference = vbc_reference_preparation(lattice,template.weights,lambda,theta)
            full = full_mpo_matrix(twisted_exchange_mpo(sites,prepared,theta),sites)[indices,indices]
            restored = full_mpo_matrix(twisted_exchange_mpo(sites,removed,theta),sites)[indices,indices]
            Q = Diagonal([2count_ones(state)-9 for state in reference.basis])
            metrics = Dict("prepared_mpo_max_error"=>maximum(abs.(full-reference.H)),
                "hermiticity_max_error"=>maximum(abs.(full-full')),
                "charge_commutator_max_error"=>maximum(abs.(Q*full-full*Q)),
                "removed_mpo_max_error"=>maximum(abs.(restored-reference_nn.H)),
                "preparation_matrix_change"=>maximum(abs.(full-reference_nn.H)))
            record["metrics"][name] = metrics
            check(name*"_matrix_contract",all(isfinite(metrics[k]) && metrics[k]<=2e-12
                for k in ("prepared_mpo_max_error","hermiticity_max_error","charge_commutator_max_error","removed_mpo_max_error")))
            check(name*"_preparation_nonzero",metrics["preparation_matrix_change"]>0.01)
            check(name*"_sites_and_winding_preserved",prepared.sites==lattice.sites &&
                [(b.i,b.j,b.wy) for b in prepared.bonds]==[(b.i,b.j,b.wy) for b in lattice.bonds])
            check(name*"_exact_removal",removed.sites==lattice.sites &&
                [(b.i,b.j,b.wy,b.Jxy,b.Jz) for b in removed.bonds]==
                [(b.i,b.j,b.wy,b.Jxy,b.Jz) for b in lattice.bonds])
            kind===:hourglass && (prep_for_solve=prepared)
            persist()
        end
        # One small trajectory suffices for the shared warm-start/checkpoint
        # interface. Both motif Hamiltonians were checked independently above.
        solver = (; Q=1,sites,seed=11,nsweeps=4,maxdim=16,cutoff=0.0,noise=0.0,
            eigsolve_tol=1e-12,eigsolve_krylovdim=20,eigsolve_maxiter=10,measure_variance=false)
        preparation = run_dmrg(prep_for_solve,0.0;solver...)
        snapshot = save_checkpoint(joinpath(output,"prepared"),preparation;
            baseline=preparation,theta_path=[0.0],status=:trial)
        saved = load_checkpoint(snapshot,prep_for_solve;Q=1,sites,status=:trial,expected_theta=0.0)
        rejected = try
            load_checkpoint(snapshot,lattice;Q=1,sites,status=:trial)
            false
        catch exception
            exception isa ArgumentError || rethrow()
            occursin("checkpoint model or geometry mismatch",sprint(showerror,exception))
        end
        check("prepared_checkpoint_rejected_as_NN",rejected)
        check("prepared_checkpoint_overlap",abs(abs(inner(saved.psi,preparation.psi))-1)<=1e-12)
        relaxed = run_dmrg(lattice,0.0;solver...,psi0=saved.psi)
        reference = reference_hamiltonian(1,3,0.0;Q=1)
        amplitudes = sector_amplitudes(relaxed.psi,sites,reference.basis)
        energy_error = abs(relaxed.energy-eigmin(Hermitian(reference.H)))/9
        residual = norm(reference.H*amplitudes-relaxed.energy*amplitudes)
        record["metrics"]["warm_start"] = Dict("energy_error_per_site"=>energy_error,
            "NN_residual_norm"=>residual,"measured_truncation_errors"=>relaxed.max_truncation_errors,
            "total_sz_error"=>abs(sum(relaxed.sz)-0.5),
            "prepared_energy"=>preparation.energy,"relaxed_NN_energy"=>relaxed.energy)
        check("released_state_matches_NN_ED",isfinite(energy_error) && energy_error<1e-8 && residual<2e-6)
        check("warm_start_complex_charge_sites",flux(relaxed.psi)==QN("Sz",1) &&
            all(eltype(t)==ComplexF64 for t in relaxed.psi) &&
            all(siteind(relaxed.psi,i)==sites[i] for i in eachindex(sites)))
        final_snapshot = save_checkpoint(joinpath(output,"released"),relaxed;
            baseline=relaxed,theta_path=[0.0],status=:trial)
        released = load_checkpoint(final_snapshot,lattice;Q=1,sites,status=:trial,expected_theta=0.0)
        check("released_checkpoint_matches_NN",abs(abs(inner(released.psi,relaxed.psi))-1)<=1e-12)
        record["status"] = "passed"
    catch exception
        record["status"] = "failed"
        record["exception_type"] = string(typeof(exception))
        rethrow()
    finally
        record["wall_elapsed_seconds"] = (time_ns()-VBC_CHECK_START)/1e9
        record["sources_unchanged"] = all(vbc_check_hash(joinpath(VBC_CHECK_ROOT,p))==h for (p,h) in hashes)
        record["sources_unchanged"] || (record["status"]="source_changed")
        persist()
    end
    return record
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==1 || error("expected NEW_OUTPUT")
    record = vbc_preparation_checks(abspath(only(ARGS)))
    println(record["status"]," in ",record["wall_elapsed_seconds"]," seconds")
    record["status"]=="passed" || exit(1)
end
