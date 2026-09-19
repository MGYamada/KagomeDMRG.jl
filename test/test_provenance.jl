@testset "Loaded implementation provenance" begin
    package_root = pkgdir(KagomeDMRG)
    identity = KagomeDMRG._execution_identity()
    worker = joinpath(@__DIR__, "provenance_worker.jl")
    mktempdir() do directory
        package_copy = joinpath(directory, "package")
        mkpath(package_copy)
        # Copy precisely the hashed sources, including the vendored backend.
        # No production file is modified by the child-process regressions.
        for (relative, _) in identity.source_sha256
            destination = joinpath(package_copy, relative)
            mkpath(dirname(destination))
            cp(joinpath(package_root, relative), destination)
        end
        output_root = joinpath(directory, "snapshots")
        checkpoint_record = joinpath(directory, "checkpoint.txt")
        isolated_depot = joinpath(directory, "depot")
        mkpath(isolated_depot)
        depot_separator = Sys.iswindows() ? ';' : ':'
        environment = ["JULIA_DEPOT_PATH" => join([isolated_depot; DEPOT_PATH], depot_separator),
                       "JULIA_NUM_THREADS" => "1", "OPENBLAS_NUM_THREADS" => "1"]
        function child(mode)
            # Exercise the stale in-memory backend directly, without rebuilding
            # a second backend in a precompile subprocess after the deliberate edit.
            cache_flags = mode == "stale_backend" ? ["--compiled-modules=existing"] : String[]
            command = `$(Base.julia_cmd()) $cache_flags --project=$package_copy --startup-file=no --threads=1 $worker $mode $output_root $checkpoint_record`
            output = IOBuffer()
            process = run(pipeline(ignorestatus(addenv(command, environment...));
                                   stdout=output, stderr=output))
            passed = success(process)
            passed || print(String(take!(output)))
            return passed
        end
        @testset "Source drift and restart in one process" begin
            @test child("compact")
        end
        @testset "Already loaded backend rejects source drift" begin
            @test child("stale_backend")
        end
        @testset "Manifest edit invalidates the precompile cache" begin
            manifest = joinpath(package_copy, identity.manifest)
            open(manifest, "a") do io
                write(io, "\n# Isolated precompile dependency regression change.\n")
            end
            @test child("new_manifest")
        end
    end
end
