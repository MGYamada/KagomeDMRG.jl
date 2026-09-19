# Source identity is captured when the module is evaluated (including precompile).
# The active environment is captured separately in __init__, since one compiled
# package may be loaded by different projects. Never relabel an in-memory result
# using source or environment files first hashed at checkpoint-save time.
const _CHECKPOINT_SOURCES = ("Project.toml", "src/KagomeDMRG.jl",
    "src/provenance.jl", "src/lattice.jl", "src/extended_model.jl", "src/model.jl", "src/dmrg.jl",
    "src/observables.jl", "src/chirality.jl", "src/checkpoint.jl", "src/schmidt.jl",
    "src/continuation.jl")

_file_sha256(path) = bytes2hex(open(sha256, path))

function _checkpoint_source_files(root)
    paths = String[_CHECKPOINT_SOURCES...]
    vendor = joinpath(root, "vendor", "NDTensors")
    isfile(joinpath(vendor, "Project.toml")) && isdir(joinpath(vendor, "src")) ||
        throw(ArgumentError("the vendored NDTensors implementation is required"))
    push!(paths, "vendor/NDTensors/Project.toml")
    for name in ("src", "ext")
        directory = joinpath(vendor, name)
        isdir(directory) || continue
        for (parent, _, files) in walkdir(directory), file in files
            push!(paths, relpath(joinpath(parent, file), root))
        end
    end
    return sort!(paths)
end

function _checkpoint_source_identity(; register_dependencies=false)
    root = dirname(@__DIR__)
    paths = _checkpoint_source_files(root)
    if register_dependencies
        # Directory dependencies also invalidate the cache when vendored source
        # files are added or removed. Environment files are not precompile inputs.
        directories = Set([root])
        for name in ("src", "ext")
            directory = joinpath(root, "vendor", "NDTensors", name)
            isdir(directory) || continue
            for (parent, _, _) in walkdir(directory)
                push!(directories, parent)
            end
        end
        for path in paths
            absolute = joinpath(root, path)
            Base.include_dependency(absolute; track_content=true)
            directory = dirname(absolute)
            while directory != root
                push!(directories, directory)
                parent = dirname(directory)
                parent != directory || error("source dependency is outside the package root")
                directory = parent
            end
        end
        for directory in sort!(collect(directories))
            Base.include_dependency(directory; track_content=true)
        end
    end
    return Tuple(path => _file_sha256(joinpath(root, path)) for path in paths)
end

function _checkpoint_environment_paths()
    project = Base.active_project()
    project !== nothing && isfile(project) ||
        throw(ArgumentError("an instantiated active project is required"))
    # Use Julia's loader selection, including version-specific names, explicit
    # manifest paths, and workspace manifests. Never fall back to this package.
    manifest = Base.project_file_manifest_path(project)
    manifest !== nothing && isfile(manifest) ||
        throw(ArgumentError("the active project needs a manifest; run Pkg.instantiate() first"))
    return (; project=abspath(project), manifest=abspath(manifest))
end

function _checkpoint_environment_identity(paths)
    # Roles prevent collisions when a custom manifest is also named Project.toml.
    hashes = ["project:" * basename(paths.project) => _file_sha256(paths.project),
              "manifest:" * basename(paths.manifest) => _file_sha256(paths.manifest)]
    return (; manifest=basename(paths.manifest), environment_sha256=Tuple(sort!(hashes)))
end

function _checkpoint_runtime_identity()
    result = Dict{String,Any}("julia" => string(VERSION), "arch" => string(Sys.ARCH),
        "kernel" => string(Sys.KERNEL), "word_size" => Sys.WORD_SIZE)
    for mod in (ITensors, ITensorMPS, ITensors.NDTensors, KrylovKit,
                LinearAlgebra, Serialization, SHA, TOML)
        result[string(nameof(mod))] = string(Base.pkgversion(mod))
    end
    return Tuple(sort!(collect(result); by=first))
end

# Keep the immutable hash tuples behind abstract fields. Propagating their
# lengths through every result NamedTuple causes excessive LLVM scalarization
# in callers that capture or merge several completed points.
struct _ExecutionIdentity
    manifest::String
    source_sha256::Tuple
    runtime::Tuple
    environment_sha256::Tuple
end

function Base.:(==)(left::_ExecutionIdentity, right::_ExecutionIdentity)
    return left.manifest == right.manifest &&
           left.source_sha256 == right.source_sha256 && left.runtime == right.runtime &&
           left.environment_sha256 == right.environment_sha256
end

const _LOADED_SOURCE_IDENTITY = _checkpoint_source_identity(; register_dependencies=true)
const _LOADED_EXECUTION_IDENTITY = Ref{_ExecutionIdentity}()
# Absolute locations are only in-memory guards, never saved as environment metadata.
const _LOADED_ENVIRONMENT_PATHS = Ref{NamedTuple{(:project, :manifest),Tuple{String,String}}}()

function _initialize_execution_identity!()
    paths = _checkpoint_environment_paths()
    environment = _checkpoint_environment_identity(paths)
    _LOADED_ENVIRONMENT_PATHS[] = paths
    _LOADED_EXECUTION_IDENTITY[] = _ExecutionIdentity(environment.manifest,
        _LOADED_SOURCE_IDENTITY, _checkpoint_runtime_identity(), environment.environment_sha256)
    _execution_identity()
    return nothing
end

function _execution_identity()
    vendor = joinpath(dirname(@__DIR__), "vendor", "NDTensors")
    loaded_backend = pkgdir(ITensors.NDTensors)
    loaded_backend !== nothing && isdir(loaded_backend) && isdir(vendor) &&
        realpath(loaded_backend) == realpath(vendor) ||
        throw(ArgumentError("the loaded NDTensors must be the package's vendored implementation"))
    source = try
        _checkpoint_source_identity()
    catch exception
        exception isa InterruptException && rethrow()
        throw(ArgumentError("loaded source is no longer readable"))
    end
    identity = _LOADED_EXECUTION_IDENTITY[]
    source == identity.source_sha256 ||
        throw(ArgumentError("source changed after KagomeDMRG was loaded; restart Julia"))
    paths, environment = try
        current = _checkpoint_environment_paths()
        (current, _checkpoint_environment_identity(current))
    catch exception
        exception isa InterruptException && rethrow()
        throw(ArgumentError("active environment is no longer readable; restart Julia"))
    end
    paths == _LOADED_ENVIRONMENT_PATHS[] &&
        environment.manifest == identity.manifest &&
        environment.environment_sha256 == identity.environment_sha256 ||
        throw(ArgumentError("active environment changed after KagomeDMRG was loaded; restart Julia"))
    isdefined(ITensors.NDTensors, :_KAGOME_LOADED_SOURCE_SHA256) ||
        throw(ArgumentError("loaded vendored NDTensors has no source identity; restart Julia"))
    vendor_prefix = joinpath("vendor", "NDTensors")
    vendor_sources = Tuple(relpath(path, vendor_prefix) => hash
        for (path, hash) in source if startswith(path, vendor_prefix))
    vendor_sources == ITensors.NDTensors._KAGOME_LOADED_SOURCE_SHA256 ||
        throw(ArgumentError("vendored NDTensors source changed after its module was loaded; restart Julia"))
    _checkpoint_runtime_identity() == identity.runtime ||
        throw(ArgumentError("runtime changed after KagomeDMRG was loaded; restart Julia"))
    return identity
end

function _require_execution_identity(point, identity=_execution_identity())
    hasproperty(point, :execution_identity) ||
        throw(ArgumentError("completed point is missing its execution identity"))
    point.execution_identity == identity ||
        throw(ArgumentError("completed point execution identity does not match the loaded implementation"))
    return identity
end

_checkpoint_runtime() = Dict(_execution_identity().runtime)

function _checkpoint_provenance(identity=_execution_identity())
    root = normpath(joinpath(@__DIR__, ".."))
    revision = try
        strip(read(Cmd(Cmd(["git", "rev-parse", "HEAD"]); dir=root), String))
    catch
        "unavailable"
    end
    dirty = try
        !isempty(read(Cmd(Cmd(["git", "status", "--porcelain"]); dir=root), String))
    catch
        "unavailable"
    end
    return Dict("source_sha256" => Dict(identity.source_sha256),
                "environment_sha256" => Dict(identity.environment_sha256),
                "manifest" => identity.manifest,
                "git_revision" => revision, "worktree_dirty" => dirty)
end
