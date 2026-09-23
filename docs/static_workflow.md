# Static run configuration and records

The public `static_run_config`, `run_static`, `load_static_run`, and
`diagnose_static_run` APIs share the configuration, solve publication, and
input-validation contract used by `examples/static_diagnostics.jl`.
The workflow prepares one zero-flux nearest-neighbor trial with configurable
`Jxy`, `Jz`, and longitudinal fields. It starts from a random complex fixed-Q
MPS, saves the completed trial, and leaves variance and other diagnostics for a
separate step. Completion does not imply convergence or branch acceptance.

Extended J1/J2/J3 models, prepared states, warm starts, and nonzero-flux paths
continue to use the existing lower-level APIs. Historical study drivers and
records retain their original protocols. No research calculation is launched
by validating or loading a configuration.

## Public workflow

Use one Julia executable and active environment for both processes. Output
parents must exist, and each output directory must be new.

```julia
using KagomeDMRG, TOML

config = static_run_config(TOML.parsefile("examples/configs/static_diagnostics_small.toml"))
solved = run_static("outputs/new-solve", config)
# solved.record: detached dictionary matching the sealed solve.toml
# solved.result: ordinary run_dmrg result, with variance === nothing
# solved.snapshot: absolute trial checkpoint directory
```

In another process using the same source and environment:

```julia
using KagomeDMRG

report = diagnose_static_run("outputs/new-solve"; output="outputs/new-diagnosis")
# Inspect report["status"], report["integrity"], and report["precision"].
```

For inspection or custom analysis, `load_static_run("outputs/new-solve")`
returns `config`, `record`, `snapshot`, `checkpoint`, and `record_sha256`.
`checkpoint` is the fully validated `load_checkpoint` result; loading performs
state and energy integrity checks, but no optimization or full diagnostics.
`diagnose_static_run` validates the envelope and calls `diagnose_checkpoint`
once; it does not load the MPS twice through `load_static_run`.

`run_static(...; progress_callback=callback)` forwards the scalar event contract
of `run_dmrg`. Callback exceptions propagate after recording their type and the
active phase. Callbacks and arbitrary exception text are not serialized.
The library records the Julia/BLAS thread counts without changing them; the
example selects one thread and enforces its own computation bounds.

## Configuration schema 1

`static_run_config` accepts a string-keyed dictionary, rejects unknown or missing
keys, and returns a detached TOML-compatible dictionary with all defaults filled.
Validating the canonical dictionary again gives the same values. Mutating an
input schedule or field vector cannot change the returned configuration.

| Key | Requirement or default |
| --- | --- |
| `schema_version` | Required integer `1` |
| `Lx`, `Ly` | Required integers; `Lx >= 1`, `Ly >= 3`; `N = 3Lx*Ly` |
| `Q` | Required integer `2Sz`, valid parity and `-N <= Q <= N` |
| `seed` | Required nonnegative machine integer |
| `nsweeps` | Required positive integer |
| `maxdim` | Required vector of positive integers; length `1:nsweeps`; last entry repeats |
| `Jxy`, `Jz` | Finite real numbers, each defaulting to `1.0` |
| `hz` | Finite real vector of length N, default all zero; Hamiltonian contains `-sum(hz .* Sz)` |
| `gauge` | String `"seam"` (default) or `"uniform"` |
| `initial_linkdim` | Positive integer, default `4` |
| `cutoff`, `noise` | Finite nonnegative reals, defaults `1e-12`, `0.0` |
| `eigsolve_tol` | Finite positive real, default `1e-12` |
| `eigsolve_krylovdim` | Integer at least 2, default `20` |
| `eigsolve_maxiter` | Positive integer, default `10` |

Boolean values are rejected in numeric fields, and floating-point values are
rejected in integer fields. Scalars are not accepted for the `maxdim` vector in
this schema. Real values are normalized to finite `Float64`, integer values to
`Int`; overflow is rejected. The workflow fixes `theta=0.0`,
`initialization="random_fixed_charge"`, and `measure_variance=false`.
These are recorded explicitly in the solve contract, not accepted as config
options. Q is always explicit; the 1/9 target uses `Q=N/9` with N divisible by 9.

The library has no small-system size, sweep, or dimension cap. The example
retains at most 18 sites, four sweeps, and dimensions 64 (including the initial
link dimension). Those bounds are development-fixture budgets, not convergence
criteria or future research budgets. The committed seven-key small config
continues to work with the materialized defaults above.

## Saved solve contract

A successful run contains:

| File | Contract |
| --- | --- |
| `config.toml` | Complete canonical configuration, including all defaults |
| `solve.toml` | Format `KagomeDMRG.static_run`, schema 1, solve status and allowlisted metadata |
| `solve.sha256` | SHA-256 of the final solve record; published only after successful solve/checkpoint publication |
| `trial/checkpoint-*/` | Existing schema-1 checkpoint (`metadata.toml`, `state.jls`, `checksums.toml`) |

`solve.toml` holds the full geometry/bond/field configuration, complete solver
settings, fixed theta and trial status, actual thread counts, solve/checkpoint/
total timings, runtime versions, and source/environment provenance. It records
the checkpoint as a relative path. Its exact `file_sha256` allowlist contains
`config.toml` and the three checkpoint files. No environment dump, hostname,
connection settings, or exception message is recorded. Provenance hash and git
revision fields retain the existing checkpoint meanings.

TOML records are replaced atomically through a same-directory temporary file.
The checkpoint uses the existing atomic directory publication. There is no
power-loss durability guarantee or transaction across the record and seal:
a process killed between their publication can leave an unsealed completed
record, which is rejected by the loader. The seal checks integrity; it is not an
authenticity signature. These are trusted local files, not an untrusted upload
or portable interchange format.

| Solve state | Meaning |
| --- | --- |
| `running` | In progress, or interrupted before an updated status was saved |
| `failed` | An exception was recorded; `active_phase`, `solver_phase`, and `exception_type` locate the failure |
| `completed_solve` with valid seal | Optimizer returned and the trial checkpoint was saved |

`diagnostic_status="not_run"` describes the deferred diagnostic step at solve
publication. It remains unchanged after later measurements. The example's
`direct-diagnostics.toml` is an optional measurement reference written after the
solve is sealed. It has its own running/failed/completed status and is not a
required solve artifact. A failure or interruption in this reference cannot
invalidate or hide the already completed solve.

## Loading, diagnosis, and compatibility

`load_static_run` checks the record format/version/status, canonical config,
reconstructed model and settings (including primitive types), exact artifact
list, hashes, runtime, source, and active environment. Artifact paths must stay
inside the run and cannot traverse symlinks. It then invokes the existing strict
checkpoint loader, which verifies checksums before deserialization, followed by
charge, site identity, norm, profile, baseline, and reconstructed energy checks.
Input hashes are checked again on return or failure.

`diagnose_static_run` shares this preflight. The output must be outside the
entire solve directory, including through symlink aliases. Invalid run envelopes
or output locations create no diagnostic output. Checkpoint-load or measurement
failures are recorded by `diagnose_checkpoint` and rethrown. Alongside
`diagnostics.toml`, `driver.toml` records the solve-record hash and whether the
pinned input files remained unchanged. Input changes cause an exception;
inspect both files when diagnosing failed work. Rerun into a new directory,
without rewriting the original solve, checkpoint, or previous diagnosis.

The new run format is distinct from the original untagged example record.
Old ad hoc records are not migrated or implicitly treated as this schema.
The checkpoint schema itself is unchanged, and strict source/environment
matching remains required. Adding workflow code changes the source identity;
older checkpoints require their original checkout/environment. Git revision
alone is not a substitute for the actual source hashes.

## Lower-level shared interfaces

These interfaces remain available for models and protocols outside this first
workflow. New drivers can consume documented result fields rather than calling
internal `_checkpoint_*` or `_execution_identity` helpers.

| Interface | Data contract used by a driver |
| --- | --- |
| `run_dmrg(lattice, theta; ...)` | Result has `psi`, `H`, `sites`, copied `lattice`/`hz`, `theta`, `gauge`, `Q`, `energy`, `local_energy`, `variance`, `sz`, per-sweep `sweep_energies`/`max_truncation_errors`, `settings`; treat attached `execution_identity` as opaque |
| `save_checkpoint(root, result; baseline, theta_path, status)` | Saves a completed point and measured zero-flux baseline; returns the new snapshot path; acceptance is caller-supplied |
| `load_checkpoint(snapshot, lattice; ...)` | Returns state/sites/model, `settings`, `diagnostics`, `metadata`, measured `baseline`, full unwrapped `theta_path`, and opaque execution identity |
| `resume_dmrg(snapshot, lattice, theta; ...)` | Starts a fresh batch from an accepted checkpoint using stored solver settings and rebuilt environments; no mid-sweep restart or branch acceptance |
| `static_diagnostics(result_or_checkpoint; ...)` | Measures without optimization; report separates measurement status, integrity, and precision |

Solver settings record initialization mode, seed, initial dimension, sweep and
maxdim schedule, cutoff, noise, eigensolver parameters, and variance choice.
`energy` is the final Hamiltonian expectation; `local_energy` and sweep values
are Ritz values. `variance === nothing` means unmeasured. Configured cutoff is
separate from measured truncation errors. Geometry, state, fields, and settings
attached to results must not be mutated and then relabeled as the old result.

API names beginning with `_` remain implementation details. Existing historical
drivers still use them; they are not all migrated by this unit. Extending the
shared workflow to prepared/continued or extended-model runs requires explicit
config and provenance contracts, not edits to historical records.

[`magnetization_curve`](magnetization_curve.md) can analyze a supplied table of
total zero-field sector energies without another solve. A shared adapter for
compatible run records is not yet implemented: callers must verify common model,
geometry, field convention, and input provenance, and preserve precision labels.

## Validation (2026-09-23)

Julia 1.13.0, one Julia/BLAS thread, fixed research dependencies:

```sh
julia --project=research --startup-file=no --threads=1 -e 'using Pkg; Pkg.test("KagomeDMRG"; allow_reresolve=false, test_args=["workflow", "provenance"])'
```

The integrated focused run passed **154/154** assertions in 159.9 seconds of
Test time; the `Pkg.test` call took 169.352 seconds, excluding outer Julia
startup. Package precompilation took about 3 seconds within that call. The
existing provenance tests also start separate processes; final per-group JIT
and numerical timings were not separated. An earlier workflow-only check passed
125 assertions in 50.1 seconds before the final failure-path additions.

The new workflow tests use one N9, one-sweep field-only solve with analytic
energy `-Q/2`, and reuse its snapshot for strict loading, successful diagnosis,
invalid-cut failure, altered records with recomputed hashes, primitive-type
changes, missing seals, artifact symlinks, forbidden output locations, and input
preservation. An early callback failure checks that arbitrary exception text is
not saved. Existing provenance checks verify source/environment drift and cached
module initialization with the new source file included. Independent review
found no blocking issue; its failure-path suggestions are covered.

The full suite, the existing fresh-process measurement equivalence fixtures, and
the refactored CLI's two-command solve/diagnose example were not rerun. Their
measurement implementation is unchanged; this unit tests the new publication and
loading contracts and moves the example's inexpensive configuration checks into
the workflow group. No current Julia 1.12 validation, research scan, N54/N72 run,
or historical result recomputation is claimed. The combined focused check still
exceeds the 60-second daily target; CI and test cost remain the next priorities.

See [static diagnostics](static_diagnostics.md) for measurement/report fields and
[test groups](../test/README.md) for validation scope and timings.
