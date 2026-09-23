#!/usr/bin/env julia
# Analyze a supplied energy table without running DMRG. See docs/magnetization_curve.md.
module MagnetizationCurveExample

using KagomeDMRG
using TOML

require(condition, message) = condition || throw(ArgumentError(message))

function require_keys(table, required, optional, label)
    require(table isa AbstractDict, "$label must be a table")
    require(all(key isa AbstractString for key in keys(table)), "$label keys must be strings")
    require(all(haskey(table, key) for key in required), "$label is missing a required key")
    allowed = Set((required..., optional...))
    require(all(key in allowed for key in keys(table)), "$label contains an unknown key")
end

"""Validate the example's allowlisted table, then construct the public report."""
function analyze_config(config)
    require_keys(config, ("schema_version", "format", "scope", "N", "fields", "sectors"),
        ("tie_atol",), "input")
    version = config["schema_version"]
    require(version isa Integer && !(version isa Bool) && version == 1,
        "schema_version must be integer 1")
    require(config["format"] == "KagomeDMRG.magnetization_energy_table", "unknown input format")
    require(config["scope"] in ("synthetic_fixture", "provided_energies"),
        "scope must be synthetic_fixture or provided_energies")
    require(config["fields"] isa AbstractVector, "fields must be an array")
    require(config["sectors"] isa AbstractVector, "sectors must be an array of tables")
    energies = Pair[]
    statuses = Dict()
    status_supplied = false
    for sector in config["sectors"]
        require_keys(sector, ("Q", "energy"), ("precision_status",), "sector")
        push!(energies, sector["Q"] => sector["energy"])
        statuses[sector["Q"]] = get(sector, "precision_status", "not_assessed")
        status_supplied |= haskey(sector, "precision_status")
    end
    curve = magnetization_curve(energies; N=config["N"], fields=config["fields"],
        tie_atol=get(config, "tie_atol", 1e-10),
        precision_status=status_supplied ? statuses : nothing)
    return Dict{String,Any}("schema_version" => 1,
        "format" => "KagomeDMRG.magnetization_example", "input_scope" => config["scope"],
        "curve" => curve)
end

# All string cells are quoted. Lists retain every candidate, including ties.
csv_cell(value::AbstractVector) = csv_cell(join(value, ";"))
csv_cell(value::AbstractString) = "\"" * replace(value, "\"" => "\"\"") * "\""
csv_cell(value) = string(value)

function write_samples(io, samples)
    columns = ("h", "minimum_energy", "status", "minimizing_Q", "near_minimizing_Q",
        "magnetization", "near_magnetization")
    println(io, join(columns, ","))
    for sample in samples
        println(io, join((csv_cell(sample[key]) for key in columns), ","))
    end
end

"""Write report.toml and samples.csv into a new output directory; never run a solver."""
function analyze(output, config_path)
    report = analyze_config(TOML.parsefile(config_path))
    destination = abspath(output)
    require(!ispath(destination) && !islink(destination), "output must be a new directory")
    require(isdir(dirname(destination)), "output parent directory must exist")
    # Render before creating the output, so invalid input/serialization leaves no files.
    report_text = sprint(io -> TOML.print(io, report; sorted=true))
    sample_text = sprint(io -> write_samples(io, report["curve"]["samples"]))
    mkdir(destination)
    write(joinpath(destination, "report.toml"), report_text)
    write(joinpath(destination, "samples.csv"), sample_text)
    return report
end

function main(args=ARGS)
    require(length(args) == 2, "usage: magnetization_curve.jl OUTPUT CONFIG")
    return analyze(args[1], args[2])
end

end # module

if abspath(PROGRAM_FILE) == @__FILE__
    MagnetizationCurveExample.main()
end
