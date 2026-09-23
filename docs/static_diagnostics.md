# Static diagnostics after checkpointing

`static_diagnostics` and `diagnose_checkpoint` share the measurement logic used
by the small development driver. They separate optimization and saving from
postprocessing. They do not run another sweep, accept a flux branch, or establish
ground-state convergence. Historical study drivers and their records retain
their original protocols; this is the first shared path for new drivers.

## Run the small driver

From the repository root, use the same Julia executable and research environment
for both commands. Each output must be new. The default config is a nine-site,
Q=1, zero-flux nearest-neighbor model with two sweeps and maxdim `[16,32]`.

```sh
julia --project=research --startup-file=no --threads=1 examples/static_diagnostics.jl solve outputs/diagnostic-solve examples/configs/static_diagnostics_small.toml
julia --project=research --startup-file=no --threads=1 examples/static_diagnostics.jl diagnose outputs/diagnostic-solve outputs/diagnostic-measurement
```

The solve directory contains an immutable trial checkpoint, `config.toml`,
`solve.toml`, its SHA-256, and `direct-diagnostics.toml`. Configuration and sealed
solve records now use the [shared static-run contract](static_workflow.md).
The direct measurement
is a validation reference for comparison with the later process. The solve
record retains `diagnostic_status = "not_run"` for the separate checkpoint
diagnostic stage; it is not rewritten after diagnosis. A failed direct-reference
measurement does not change the already completed solve or delete its saved trial.
This optional reference is separate from the sealed solve artifacts.

The diagnostic directory contains `diagnostics.toml` and `driver.toml`. The
shared workflow checks its source, environment, solve record, config, and all checkpoint file hashes
before and after measuring. Rerun into another new directory; the original solve
record, parent checkpoint, and previous diagnostic outputs remain unchanged.
Preflight errors (for example an altered solve record) abort before measurement.

The committed example config supplies `schema_version`, `Lx`, `Ly`, `Q`, `seed`,
`nsweeps`, and `maxdim`; the shared schema materializes model and solver defaults
and accepts the optional fields documented in the [configuration contract](static_workflow.md#configuration-schema-1).
The example deliberately bounds this development fixture to at
most 18 sites, four sweeps, and dimension 64 (including the initial link). These are example limits, not limits
of the library API or future research budgets. The example uses one Julia/BLAS
thread. It is not a convergence study.

## Public API

```julia
direct = static_diagnostics(result; cuts=1:result.lattice.Lx-1,
                            chirality=false)
saved = load_checkpoint(snapshot, lattice; Q, gauge, hz, status=:trial)
in_memory = static_diagnostics(saved)
report = diagnose_checkpoint(snapshot, lattice;
    output="outputs/new-diagnostic", Q, gauge, hz,
    expected_theta=0.0, expected_settings=result.settings, status=:trial)
```

`static_diagnostics` accepts a completed `run_dmrg` result or a strictly loaded
checkpoint. It validates the execution identity and copies the MPS before
measurement. `diagnose_checkpoint` calls the existing strict loader and writes
only into a newly created output directory whose parent already exists. Output
inside the checkpoint is rejected. All three checkpoint files are hashed before
and after, including on failure. Existing source/environment/runtime, geometry,
charge, field, site identity, and checkpoint checksum checks remain in force.

Both APIs use the actual oriented bonds, total Q, unwrapped theta, gauge and
longitudinal fields. Bond energies include the complex exchange phases; the
sum minus `dot(hz, sz)` is compared with total energy. The Hamiltonian is rebuilt
from these inputs. Stored energy and Sz are remeasured for consistency, then
retained as the reference convention for variance and profile comparisons.
The raw variance is `real(<H†H>) - E²`; it is not an energy error bar. A small
negative roundoff value is retained and checked, not clipped to zero.

`cuts` contains unique integer cell boundaries in `1:Lx-1`; an empty list is
valid, including the default for Lx=1. The corresponding MPS bond is `3Ly*cut`.
`chirality=true` adds the existing dressed triangle observable. No central-window
or CSL-specific geometry is built in. Correlation matrices still require O(N²)
output storage; selected-pair performance work remains future work.

`progress_callback(event)` on `static_diagnostics` receives scalar phase/status
and elapsed-time notifications. Exceptions propagate. `diagnose_checkpoint`
uses this internally to persist its current phase; it records a failure before
rethrowing load or measurement exceptions. Exception types, not arbitrary error
messages, are saved. A forcibly terminated process may leave `status="running"`;
that status is not evidence of completion.

## Report contract (schema 1)

Reports are TOML-compatible dictionaries with
`format = "KagomeDMRG.static_diagnostics"`. They contain selected model/configuration,
solver settings, source/environment provenance, runtime versions, and per-phase
timing. They never serialize callbacks or all environment variables.

| Field | Meaning |
| --- | --- |
| `status` | `running`, `completed`, or `failed`; only the wrapper persists intermediate/failure records |
| `integrity.status` | `not_evaluated`, `passed`, or `failed`; consistency with the actual model and state |
| `precision.status` | `not_evaluated`, `incomplete`, `unmet`, or `passed`; separate from execution and integrity |
| `measurements` | Energy, Sz/column profiles, bond energy, real/imaginary correlation rows, raw H†H/variance, Schmidt data, original sweep/truncation history, optional chirality |
| `measurement_seconds` | Timings for Hamiltonian construction, state checks, correlations, variance, Schmidt and optional chirality, including first-call compilation where applicable |
| `checkpoint` / `checkpoint_unchanged` | Wrapper-only input identity and final byte-preservation check |

`completed` can coexist with failed integrity or unmet precision. Load/measurement
failure has `status="failed"`, not a successful empty set of diagnostics.
Measurements preserve `settings.measure_variance=false` when variance was deferred
at solve time; computing it later does not relabel the original solver run.

The five [static precision criteria](research/p4_static_plateau_strategy.md#精度と解釈)
are unchanged: energy/site 1e-6, Sz change 1e-4, bond change 1e-4,
absolute variance/site 1e-5, and final measured truncation error 1e-6.
Pass a previous same-model/theta completed, integrity-passed report as `reference`
to either API to evaluate parent energy, Sz and bond changes. Without it those
three criteria appear in `precision.missing`. Any failed measured criterion gives
`unmet`; otherwise missing criteria give `incomplete`. Failed integrity prevents
precision certification (`not_evaluated`).
References must have the same source/environment/runtime identity as the point.
The output binds the entire reference report by SHA-256 of its sorted TOML
representation (`reference.sha256`); this includes its recorded settings,
measurements, provenance and timings. With no reference, `reference.provided`
is false. A historical report is not silently adopted into a current comparison.

With a reference, energy change is the maximum of parent difference, final
solver/expectation difference, and the final two-sweep difference when available.
A one-sweep call uses the first two without inventing an internal difference;
`within_batch_sweep_change_measured` records that distinction. Precision across
different chi schedules is labeled `chi_change`. `stationarity_passed` requires
a fixed common chi and passing all five conditions. Even that does not establish
the ground state or a phase. Repeating measurements of the same state alone is
not a convergence study.

## Validation

The `diagnostics` test group checks direct versus a fresh Julia process using the
same active environment, immutable inputs on success/failure/rerun, schema and
status semantics, and independent spin-basis measurements of a complex state at
nonzero flux with nondefault Q, uniform gauge and longitudinal fields. An analytic
two-column state checks geometric Schmidt cuts. These are development checks;
N54/N72 research calculations are not part of this workflow.

```sh
julia --project=research --startup-file=no --threads=1 -e 'using Pkg; Pkg.test("KagomeDMRG"; test_args=["diagnostics"], allow_reresolve=false)'
```

The checkpoint format remains a strict local restart format. Adding or changing
measurement source changes its source identity. Old checkpoints must still be
loaded with their original source and dependency environment; no migration or
provenance bypass is provided here.

### Small-driver check (2026-09-23)

Before extracting the shared run contract, the N9 example was run in separate Julia 1.13.0 processes with one
Julia/BLAS thread, using the fixed research environment. Direct and saved-state
measurement values matched exactly in this run (maximum absolute difference
0.0); the solve record and all pinned artifacts retained their hashes. The
diagnostic result was `completed`, integrity `passed`, precision `incomplete`
because no prior state was supplied. The raw variance/site was
`-1.578983857244667e-15`, within the existing roundoff allowance.

The timed driver calls took 45.373 seconds for solve/save/direct reference and
30.615 seconds for separate diagnostics. Julia reported 98.67% and 99.39%
compilation time respectively; these call timings exclude outer process startup
and are not steady-state performance benchmarks. Local artifacts are retained
under `outputs/codebase-diagnostics-20260923/`. No research-size state was run.

That initial diagnostics implementation's integrated package suite passed **1,855/1,855** assertions in
300.3 seconds of Test time (`Pkg.test` call: 305.561 seconds, excluding outer
Julia startup). It includes all three saved-state fixtures in one child process
and the existing checkpoint/provenance regressions. Independent code review
checked equations, statuses and input preservation. The approximately five-minute
suite still exceeds the daily three-minute target; see the
[test guide](../test/README.md) for the timing limitation. This change was tested
on Julia 1.13.0; no new Julia 1.12 test result is claimed.
These timings and artifacts describe that earlier implementation, not a rerun of
the refactored driver. Current run-contract checks are recorded in the
[workflow guide](static_workflow.md) and [test guide](../test/README.md).
