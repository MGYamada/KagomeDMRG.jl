# Run only from test_provenance.jl, in a disposable package copy and process.
using Test
using LinearAlgebra
using ITensors
using ITensorMPS
include("checkpoint_helpers.jl")

BLAS.set_num_threads(1)
ITensors.disable_threaded_blocksparse()

if first(ARGS) == "stale_backend"
    # Load the backend first, then change its implementation before KagomeDMRG
    # evaluates its own source identity. A matching backend path is insufficient.
    source = joinpath(pkgdir(ITensors.NDTensors), "src", "blocksparse", "truncate_spectrum.jl")
    original = read(source, String)
    try
        write(source, original * "\n# Isolated already-loaded backend regression.\n")
        @testset "Backend loaded before KagomeDMRG source validation" begin
            _checkpoint_test_rejection("vendored NDTensors source changed";
                                       exception_type=Exception) do
                @eval using KagomeDMRG
            end
        end
    finally
        write(source, original)
    end
    exit()
end

using KagomeDMRG

function provenance_worker_edit(action, path)
    original = read(path, String)
    try
        write(path, original * "\n# Isolated source-identity regression change.\n")
        action()
    finally
        write(path, original)
    end
end

mode, output_root, checkpoint_record = ARGS
source = joinpath(pkgdir(KagomeDMRG), "src", "dmrg.jl")

@testset "Checkpoint execution identity worker: $mode" begin
    if mode == "compact"
        # Capture must happen at module load, even before the first optimizer.
        provenance_worker_edit(source) do
            _checkpoint_test_rejection("changed after KagomeDMRG was loaded") do
                _checkpoint_fixed_field_result()
            end
        end
        result = _checkpoint_fixed_field_result()
        checkpoint = save_checkpoint(output_root, result; baseline=result,
            theta_path=[0.0], status=:accepted)
        write(checkpoint_record, checkpoint)
        loaded = load_checkpoint(checkpoint, result.lattice; hz=result.hz)
        @test loaded.execution_identity == result.execution_identity
        @test loaded.baseline.execution_identity == result.execution_identity
        @test Dict(result.execution_identity.source_sha256) ==
            loaded.metadata["provenance"]["source_sha256"]
        @test Dict(result.execution_identity.environment_sha256) ==
            loaded.metadata["provenance"]["environment_sha256"]
        @test !isfile(joinpath(pkgdir(KagomeDMRG), "Manifest.toml"))
        @test !isfile(joinpath(pkgdir(KagomeDMRG), "Manifest-v1.13.toml"))
        # A mutable environment is checked independently of immutable loaded
        # source. Reuse the completed point without another numerical solve.
        paths = KagomeDMRG._checkpoint_environment_paths()
        for environment_file in (paths.project, paths.manifest)
            provenance_worker_edit(environment_file) do
                for action in (
                    () -> _checkpoint_fixed_field_result(),
                    () -> save_checkpoint(output_root, result; baseline=result,
                        theta_path=[0.0], status=:accepted),
                    () -> load_checkpoint(checkpoint, result.lattice; hz=result.hz),
                    () -> resume_dmrg(checkpoint, result.lattice, 0.1; hz=result.hz),
                    () -> continue_flux(result.lattice, [0.1]; start=result,
                        output_root, policy=_checkpoint_test_policy(),
                        initial_step=0.1, min_step=0.01, hz=result.hz))
                    _checkpoint_test_rejection(action, "active environment changed")
                end
            end
            @test KagomeDMRG._execution_identity() == result.execution_identity
        end
        # Equal bytes at a different active location must not relabel already
        # loaded code. Absolute paths remain private and never enter metadata.
        mktempdir() do alternate
            alternate_project = joinpath(alternate, basename(paths.project))
            cp(paths.project, alternate_project)
            cp(paths.manifest, joinpath(alternate, basename(paths.manifest)))
            try
                Base.set_active_project(alternate_project)
                _checkpoint_test_rejection("active environment changed") do
                    KagomeDMRG._execution_identity()
                end
            finally
                Base.set_active_project(paths.project)
            end
        end
        @test KagomeDMRG._execution_identity() == result.execution_identity
        provenance_worker_edit(source) do
            for (label, action) in (
                ("save", () -> save_checkpoint(output_root, result; baseline=result,
                    theta_path=[0.0], status=:accepted)),
                ("load", () -> load_checkpoint(checkpoint, result.lattice; hz=result.hz)),
                ("resume", () -> resume_dmrg(checkpoint, result.lattice, 0.1; hz=result.hz)),
                ("continue", () -> continue_flux(result.lattice, [0.1]; start=result,
                    output_root, policy=_checkpoint_test_policy(),
                    initial_step=0.1, min_step=0.01, hz=result.hz)))
                @testset "$label rejects changed source" begin
                    _checkpoint_test_rejection(action, "changed after KagomeDMRG was loaded")
                end
            end
        end
        @test load_checkpoint(checkpoint, result.lattice; hz=result.hz).theta == 0.0
        # Preserve identity through ordinary restart and continuation restart
        # in the same worker, sharing the first solver compilation.
        (; lattice, hz) = _checkpoint_fixed_field_model()
        checkpoint = read(checkpoint_record, String)
        loaded = load_checkpoint(checkpoint, lattice; hz)
        resumed = resume_dmrg(checkpoint, lattice, 0.1; hz)
        @test resumed.execution_identity == loaded.execution_identity
        resaved = save_checkpoint(output_root, resumed; baseline=loaded.baseline,
            theta_path=[0.0, 0.1], status=:accepted)
        @test load_checkpoint(resaved, lattice; hz).baseline.execution_identity ==
            loaded.baseline.execution_identity
        trajectory = continue_flux(lattice, [0.1]; start=checkpoint,
            output_root, policy=_checkpoint_test_policy(), initial_step=0.1,
            min_step=0.01, hz)
        @test trajectory.status === :completed
        @test load_checkpoint(trajectory.last_checkpoint, lattice; hz).execution_identity ==
            loaded.execution_identity
    elseif mode == "new_manifest"
        # The parent changed the environment after the package was compiled.
        # __init__ must capture the current files even with a reusable cache.
        result = _checkpoint_fixed_field_result()
        identity = result.execution_identity
        paths = KagomeDMRG._checkpoint_environment_paths()
        @test Dict(identity.environment_sha256)["manifest:" * identity.manifest] ==
            bytes2hex(sha256(read(paths.manifest)))
        @test Dict(identity.environment_sha256)["project:" * basename(paths.project)] ==
            bytes2hex(sha256(read(paths.project)))
        _checkpoint_test_rejection("checkpoint source or dependency manifest mismatch") do
            load_checkpoint(read(checkpoint_record, String), result.lattice; hz=result.hz)
        end
        checkpoint = save_checkpoint(output_root, result; baseline=result,
            theta_path=[0.0], status=:accepted)
        @test load_checkpoint(checkpoint, result.lattice; hz=result.hz).execution_identity == identity
    else
        error("unknown provenance worker mode")
    end
end
