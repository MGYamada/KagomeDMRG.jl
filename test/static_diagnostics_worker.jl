# A fresh process must load the same active environment as its parent test.
using LinearAlgebra
using ITensors
using KagomeDMRG

BLAS.set_num_threads(1)
ITensors.disable_threaded_blocksparse()
length(ARGS) > 0 && length(ARGS) % 3 == 0 ||
    error("expected label/checkpoint/output triples")
for offset in 1:3:length(ARGS)
    label, snapshot, output = ARGS[offset:offset+2]
    if label == "optimized"
        diagnose_checkpoint(snapshot, kagome_cylinder(1, 3); output,
            status=:trial, expected_theta=0.0, Q=1, chirality=true)
    elseif label == "complex"
        diagnose_checkpoint(snapshot, kagome_cylinder(1, 3; Jxy=0.7, Jz=1.3); output,
            status=:trial, expected_theta=0.37, Q=-1, gauge=:uniform,
            hz=collect(range(-0.2, 0.3; length=9)))
    elseif label == "schmidt"
        diagnose_checkpoint(snapshot, kagome_cylinder(2, 3; Jxy=0, Jz=0); output,
            status=:trial, expected_theta=0.0, Q=0,
            hz=[fill(1.0, 9); fill(-1.0, 9)])
    else
        error("unknown static diagnostic fixture")
    end
end
