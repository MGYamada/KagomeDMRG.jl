using SHA
using Serialization
using TOML

function _checkpoint_test_file_hashes(directory)
    return Dict(relpath(joinpath(root, name), directory) =>
                bytes2hex(sha256(read(joinpath(root, name))))
                for (root, _, files) in walkdir(directory) for name in files)
end

# Recompute the public integrity envelope after an intentional data change.
# This distinguishes semantic validation from accidental corruption detection.
function _checkpoint_test_rehash!(directory, filename)
    path = joinpath(directory, filename)
    checksum_path = joinpath(directory, "checksums.toml")
    checksums = TOML.parsefile(checksum_path)
    checksums["files"][filename]["sha256"] = bytes2hex(sha256(read(path)))
    checksums["files"][filename]["bytes"] = filesize(path)
    open(checksum_path, "w") do io
        TOML.print(io, checksums; sorted=true)
    end
    return directory
end

function _checkpoint_test_edit_metadata!(edit, directory)
    path = joinpath(directory, "metadata.toml")
    metadata = TOML.parsefile(path)
    edit(metadata)
    open(path, "w") do io
        TOML.print(io, metadata; sorted=true)
    end
    return _checkpoint_test_rehash!(directory, "metadata.toml")
end

function _checkpoint_test_rejection(action, expected_message; exception_type=ArgumentError)
    exception = try
        action()
        nothing
    catch error
        error
    end
    @test exception isa exception_type
    @test exception isa Exception && occursin(expected_message, sprint(showerror, exception))
end

# A deliberately altered model with an analytic, flux-independent product
# ground state. It is a protocol control, not a Heisenberg phase prediction.
function _checkpoint_fixed_field_model(Lx=1; Q=nothing)
    lattice = kagome_cylinder(Lx, 3; Jxy=0.0, Jz=0.0)
    N = nsites(lattice)
    sector = target_sector(N; Q)
    nup, ndown = sector.Nup, sector.Ndown
    labels = [fill("Up", nup); fill("Dn", ndown)]
    hz = [fill(1.0, nup); fill(-1.0, ndown)]
    return (; lattice, hz, labels)
end

function _checkpoint_fixed_field_result(Lx=1; Q=nothing)
    (; lattice, hz, labels) = _checkpoint_fixed_field_model(Lx; Q)
    sites = spin_sites(lattice)
    psi0 = MPS(ComplexF64, sites, labels)
    return run_dmrg(lattice, 0.0; Q, sites, psi0, hz, nsweeps=2, maxdim=2,
                    cutoff=0.0, noise=0.0)
end

function _checkpoint_test_policy(; kwargs...)
    thresholds = (; min_overlap=0.0, max_density_change=1.0,
        max_entropy_change=10.0, max_schmidt_change=9.0,
        max_variance=1e-9, max_truncation_error=1e-10,
        max_sweep_energy_change=1e-9, max_cut_spread=1e-9,
        consistency_tol=1e-9)
    return FluxPolicy(; merge(thresholds, (; kwargs...))...)
end

# Reuse the saved product eigenstate only in control-flow tests. Checkpoint
# validation still independently rebuilds H and checks its energy and profile.
# Copy the *supplied* state, so corruption exposes a missing rollback/reload.
function _checkpoint_fixed_field_step(saved, lattice, theta; gauge, hz, outputlevel)
    all(b -> iszero(b.Jxy) && iszero(b.Jz), lattice.bonds) ||
        error("the control solver requires zero exchange couplings")
    hz == saved.hz && gauge == saved.gauge || error("fixed-field control changed")
    return merge(KagomeDMRG._flux_result(saved),
                 (; theta=Float64(theta), psi=deepcopy(saved.psi)))
end
