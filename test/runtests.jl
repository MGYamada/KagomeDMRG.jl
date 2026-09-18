using Test
using LinearAlgebra
using KagomeDMRG
using ITensors
using ITensorMPS

BLAS.set_num_threads(1)
ITensors.disable_threaded_blocksparse()

include("reference_ed.jl")
include("itensor_helpers.jl")
include("test_lattice.jl")
include("test_reference_ed.jl")
include("test_itensor.jl")
include("test_observables.jl")
