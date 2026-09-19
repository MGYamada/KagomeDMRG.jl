@testset "Loaded implementation provenance" begin
    package_root = pkgdir(KagomeDMRG)
    identity = KagomeDMRG._execution_identity()
    worker = joinpath(@__DIR__, "provenance_worker.jl")
    mktempdir() do directory
        # A valid explicit manifest location can share the project's basename.
        # Both roles must survive conversion to saved metadata dictionaries.
        custom = joinpath(directory, "custom")
        mkpath(joinpath(custom, "locks"))
        project = joinpath(custom, "Project.toml")
        manifest = joinpath(custom, "locks", "Project.toml")
        write(project, "manifest = \"locks/Project.toml\"\n")
        write(manifest, "manifest_format = \"2.0\"\n")
        environment = KagomeDMRG._checkpoint_environment_identity((; project, manifest))
        hashes = Dict(environment.environment_sha256)
        @test length(hashes) == 2
        @test hashes["project:Project.toml"] == bytes2hex(sha256(read(project)))
        @test hashes["manifest:Project.toml"] == bytes2hex(sha256(read(manifest)))
        package_copy = joinpath(directory, "package")
        mkpath(package_copy)
        # Copy precisely the loaded package sources, including the vendored
        # backend, then add a separate reproducible execution environment.
        # Deliberately leave the package root without any Manifest.
        for (relative, _) in identity.source_sha256
            destination = joinpath(package_copy, relative)
            mkpath(dirname(destination))
            cp(joinpath(package_root, relative), destination)
        end
        environment_copy = joinpath(package_copy, "research")
        mkpath(environment_copy)
        for name in ("Project.toml", "Manifest.toml", "Manifest-v1.13.toml")
            cp(joinpath(package_root, "research", name), joinpath(environment_copy, name))
        end
        @test !isfile(joinpath(package_copy, "Manifest.toml"))
        @test !isfile(joinpath(package_copy, "Manifest-v1.13.toml"))
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
            cache_flags = mode == "stale_backend" ? ["--compiled-modules=existing"] :
                          mode == "new_manifest" ? ["--compiled-modules=strict"] : String[]
            command = `$(Base.julia_cmd()) $cache_flags --project=$environment_copy --startup-file=no --threads=1 $worker $mode $output_root $checkpoint_record`
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
        @testset "Fresh process captures the active environment after precompile" begin
            manifest = Base.project_file_manifest_path(joinpath(environment_copy, "Project.toml"))
            open(manifest, "a") do io
                write(io, "\n# Isolated precompile dependency regression change.\n")
            end
            @test child("new_manifest")
        end
    end
end
