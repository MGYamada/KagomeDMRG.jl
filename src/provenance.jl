# This identity is evaluated with the module, including when it is precompiled.
# Never replace it with hashes first collected at save time: files on disk may
# already differ from the implementation that produced an in-memory result.
const _CHECKPOINT_SOURCES = ("Project.toml", "src/KagomeDMRG.jl",
    "src/provenance.jl", "src/lattice.jl", "src/extended_model.jl", "src/model.jl", "src/dmrg.jl",
    "src/observables.jl", "src/checkpoint.jl", "src/schmidt.jl",
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
    versioned = "Manifest-v$(VERSION.major).$(VERSION.minor).toml"
    manifest = isfile(joinpath(root, versioned)) ? versioned : "Manifest.toml"
    paths = [_checkpoint_source_files(root); manifest]
    if register_dependencies
        # Directory dependencies also invalidate a cache when a version-specific
        # manifest or a new vendored source file is added or removed.
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
    return (; manifest, source_sha256=Tuple(path => _file_sha256(joinpath(root, path))
                                           for path in paths))
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
end

function Base.:(==)(left::_ExecutionIdentity, right::_ExecutionIdentity)
    return left.manifest == right.manifest &&
           left.source_sha256 == right.source_sha256 && left.runtime == right.runtime
end

const _LOADED_EXECUTION_IDENTITY = let
    source = _checkpoint_source_identity(; register_dependencies=true)
    _ExecutionIdentity(source.manifest, source.source_sha256, _checkpoint_runtime_identity())
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
        throw(ArgumentError("loaded source or dependency manifest is no longer readable"))
    end
    source.manifest == _LOADED_EXECUTION_IDENTITY.manifest &&
        source.source_sha256 == _LOADED_EXECUTION_IDENTITY.source_sha256 ||
        throw(ArgumentError("source or dependency manifest changed after KagomeDMRG was loaded; restart Julia"))
    isdefined(ITensors.NDTensors, :_KAGOME_LOADED_SOURCE_SHA256) ||
        throw(ArgumentError("loaded vendored NDTensors has no source identity; restart Julia"))
    vendor_prefix = joinpath("vendor", "NDTensors")
    vendor_sources = Tuple(relpath(path, vendor_prefix) => hash
        for (path, hash) in source.source_sha256 if startswith(path, vendor_prefix))
    vendor_sources == ITensors.NDTensors._KAGOME_LOADED_SOURCE_SHA256 ||
        throw(ArgumentError("vendored NDTensors source changed after its module was loaded; restart Julia"))
    _checkpoint_runtime_identity() == _LOADED_EXECUTION_IDENTITY.runtime ||
        throw(ArgumentError("runtime changed after KagomeDMRG was loaded; restart Julia"))
    return _LOADED_EXECUTION_IDENTITY
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
                "manifest" => identity.manifest,
                "git_revision" => revision, "worktree_dirty" => dirty)
end
