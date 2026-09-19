const TEST_GROUPS = [
    "lattice" => ["test_lattice.jl", "test_extended_model.jl"],
    "ed" => ["test_reference_ed.jl", "test_reference_eigensolve.jl", "test_extended_reference.jl"],
    "dmrg" => ["test_itensor.jl"],
    "observables" => ["test_observables.jl", "test_schmidt.jl", "test_chirality.jl"],
    "checkpoint" => ["test_checkpoint.jl"],
    "truncation" => [joinpath("..", "vendor", "NDTensors", "test",
                             "test_truncation_selection.jl"),
                     "test_truncation_guard.jl", "test_qn_calibration.jl"],
    "continuation" => ["test_continuation.jl"],
    "provenance" => ["test_provenance.jl"],
]
const SELECTED_GROUPS = ARGS
if SELECTED_GROUPS == ["--help"]
    println("Usage: test/runtests.jl [group ...]")
    println("Groups: ", join(first.(TEST_GROUPS), ", "))
    println("No groups selects all regression tests.")
    exit()
end
unknown = setdiff(SELECTED_GROUPS, first.(TEST_GROUPS))
isempty(unknown) || error("Unknown test group(s): $(join(unknown, ", ")). Use --help.")

using Test
using LinearAlgebra
using KagomeDMRG
using ITensors
using ITensorMPS

BLAS.set_num_threads(1)
ITensors.disable_threaded_blocksparse()

include("reference_ed.jl")
include("reference_extended_ed.jl")
include("reference_eigensolve.jl")
include("itensor_helpers.jl")
include("checkpoint_helpers.jl")
include("truncation_helpers.jl")

@testset "KagomeDMRG" begin
    for (group, files) in TEST_GROUPS
        isempty(SELECTED_GROUPS) || group in SELECTED_GROUPS || continue
        @info "Testing $group"
        @testset "$group" begin
            for file in files
                include(file)
            end
        end
    end
end
