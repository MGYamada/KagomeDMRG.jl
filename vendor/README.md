# Pinned NDTensors correction

`NDTensors/` is a local copy of the MIT-licensed NDTensors 0.4.31 package
from the previously locked registry tree
`311282faab99a123706a738fcceb9caba81225ae`. Its original license is retained.
The local version is `0.4.31+1`; both Julia manifests and the root project's
`[sources]` entry select this directory. It is a local research correction,
not an upstream release. No files in the user's shared Julia package cache
are modified.

The numerical correction is limited to block-sparse SVD/Hermitian-eigen rank allocation in
`src/blocksparse/linearalgebra.jl` and its new helper
`src/blocksparse/truncate_spectrum.jl`. The dense decompositions, DMRG sweeps,
Hamiltonian, charge conventions, and noise construction remain upstream.
The general non-Hermitian eigen path also remains upstream: its within-block
eigenvalues do not satisfy the descending-weight premise of this selector.

Previously the global truncated spectrum selected a rank, but block allocation
reconstructed that choice from a scalar threshold. Equal and nearly equal
weights could therefore produce a different retained rank and loss. The local
correction carries the selected global entries into per-block counts directly.
For the default `min_blockdim=0`, the existing global cutoff/rank decision and
error normalization are retained. Exact ties are ordered by block coordinate
and then by within-block index, and `maxdim` is a hard cap. This can split a
degenerate subspace; it does not preserve an entire degeneracy multiplet.

With a positive `min_blockdim`, required per-block prefixes consume the same
global budget. Remaining slots are assigned by descending weight. Impossible
per-block minima exceeding `maxdim` raise `ArgumentError`. The spectrum and
discarded weight describe the states actually retained, including these
constraints. This explicitly replaces the old behavior that could exceed the
global cap without updating the reported spectrum.

The package's complete `src/`, `ext/`, and `Project.toml` contents are part of
KagomeDMRG's execution identity. A separate immutable fingerprint captured
inside NDTensors at module load detects edits made after the backend was loaded
but before KagomeDMRG was loaded. This requires only the SHA standard library;
it is separate from the numerical correction. Loading also checks that this is
the backend actually in use. Work from the repository's project environment. An application
using KagomeDMRG as a dependency must also select this same local NDTensors
source; the registry backend is rejected. A later upstream release must pass
the independent calibration before replacing this dependency.

Primary implementation and API references:

- [NDTensors block decompositions](https://github.com/ITensor/ITensors.jl/blob/main/NDTensors/src/blocksparse/linearalgebra.jl)
- [ITensor factorization semantics](https://docs.itensor.org/ITensors/stable/ITensorType.html)

See `test/test_qn_calibration.jl` for the independent dense spin-basis tests,
and `examples/validate_finite_chi.jl` for replay of a captured failure and
the bounded 18-site comparison. CPU ComplexF64 and Hermitian density matrices
are the research validation scope; no GPU performance or generic non-Hermitian
factorization claim is made.
