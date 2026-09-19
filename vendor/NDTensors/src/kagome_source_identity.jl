# KagomeDMRG local provenance support. Keep this tuple immutable so a process
# which loaded NDTensors before editing its sources cannot relabel that backend
# when it later loads KagomeDMRG.
import SHA

const _KAGOME_LOADED_SOURCE_SHA256 = let
    root = dirname(@__DIR__)
    paths = ["Project.toml"]
    directories = String[root]
    for name in ("src", "ext")
        directory = joinpath(root, name)
        isdir(directory) || continue
        for (parent, _, files) in walkdir(directory)
            push!(directories, parent)
            append!(paths, relpath(joinpath(parent, file), root) for file in files)
        end
    end
    sort!(paths)
    for relative in paths
        Base.include_dependency(joinpath(root, relative); track_content=true)
    end
    for directory in directories
        Base.include_dependency(directory; track_content=true)
    end
    Tuple(relative => bytes2hex(open(SHA.sha256, joinpath(root, relative)))
          for relative in paths)
end
