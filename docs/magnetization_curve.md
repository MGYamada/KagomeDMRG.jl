# Magnetization curves from sector energies

`magnetization_curve` constructs a zero-temperature finite-system magnetization
curve from supplied total zero-field energies. It is an analysis API: it does
not run DMRG, load checkpoints, or assess convergence. A versioned TOML example
exports the report and a CSV table without adding plotting dependencies.

## Conventions and input

For N spin-1/2 sites, integer tensor charge is `Q = 2Sz`, with allowed charges
`-N:2:N`. Physical magnetization is `M = Q/2`, saturation is `Msat = N/2`, and
the normalized magnetization is `M/Msat = Q/N`. The uniform-field Hamiltonian
has the shift `-h*Q/2`, so the curve minimizes

```text
E_h(Q) = E0(Q) - h*Q/2
```

Each energy must be a **total** energy of the same model, geometry, couplings,
and fixed background fields, with the uniform field h being scanned removed.
Do not supply energy per site or combine unrelated runs. At fixed Q, scanning h
needs no additional DMRG solve. The API accepts positive integer N with N+1
representable as `Int`, including small
artificial fixtures; the kagome 1/9 target separately requires N divisible by 9
and `Q=N/9`.

```julia
using KagomeDMRG

# Artificial example, not a kagome calculation.
energies = [-4 => 4.0, -2 => 1.0, 0 => 0.0, 2 => 1.0, 4 => 4.0]
curve = magnetization_curve(energies; N=4, fields=[0.0, 1.0, 2.0])
curve["samples"][2]["minimizing_Q"] # [0, 2]: preserve both at h=1
```

The input is a dictionary or iterable of distinct `Q=>E0` pairs. Charges must be
integers with the correct parity and range; at least one is required. Fields
are a vector, range, or tuple and keep their requested order, including repeats.
Numeric energies, fields, and `tie_atol` are normalized to finite `Float64`;
Booleans, nonfinite values, and overflow are rejected. N and Q must be integers,
not floating-point approximations or Booleans. `tie_atol` defaults to `1e-10`
and must be nonnegative.

Every supplied sector competes against **all** other supplied sectors. This
handles jumps over intermediate Q values and multi-sector crossings; an adjacent
sector estimate alone need not give the correct interval. Negative sectors are
not inferred from positive ones, and no spin-reversal symmetry is assumed.

## Report schema 1

The returned TOML-compatible dictionary has format
`KagomeDMRG.magnetization_curve`, `schema_version=1`, and
`scope="provided_sector_energies"`.

| Field | Meaning |
| --- | --- |
| `sectors` | Sorted Q, normalized total `energy`, physical `M`, normalized `magnetization`, and `precision_status` |
| `coverage` | `compared_Q`, every `missing_Q`, `complete`, and `spin_reversal_assumed=false` |
| `intervals` | One entry per supplied Q, with `h_lower`, `h_upper`, and `status` |
| `transitions` | Ordered crossing fields `h`, `field_is_exact`, `from_Q`, `to_Q`, all `coexisting_Q` and their normalized magnetizations |
| `samples` | Requested `h`, minimum field-dependent total energy, strict and nearby Q lists, magnetizations, and sample `status` |

The report also records N, charge/energy/field conventions, `tie_atol`,
`precision_status_source`, and `numerical_quality="not_assessed"`.

An interval marked `stable` has positive open width; its finite endpoints are
included as degenerate minima at exact crossings. `point` means the sector is
a minimum at a single field but has no open interval. `skipped` means it is
never a minimum; its lower constraint exceeds its upper constraint, so these
bounds do not describe a physical stable interval. Endpoints can be `-Inf` or
`Inf` when comparison with the supplied competitors is unbounded. A single
supplied sector is therefore stable for all h within this restricted comparison.

`coverage.complete=true` means every allowed Q was supplied. It says nothing
about whether those energies are converged ground-state energies. Missing
sectors can change or eliminate apparently stable intervals, including those
within the interior of the provided Q range. A partial table never establishes
the full magnetization curve, even over positive h alone.

## Arithmetic, crossings, and near ties

After normalization, the algorithm compares exact rational representations of
the `Float64` inputs. This preserves the envelope topology for those represented
numbers, without using a tolerance to erase a narrow but positive interval. It
does **not** certify that input numerical energies equal physical ground-state
energies, and digits beyond Float64 precision are not retained from input types.

Finite interval bounds, transition fields, and sampled minimum energies are
returned as Float64. Nonfinite conversion or distinct bounds/transitions that
collapse to the same Float64 value is rejected. Infinity denotes a genuinely
unbounded endpoint, not arithmetic overflow.

For a transition, `coexisting_Q` refers to the exact rational crossing and
includes sectors stable only at that point. `from_Q` and `to_Q` are its smallest
and largest coexisting charges. `field_is_exact=false` means the returned h was
rounded to Float64. Sampling that rounded h can give a unique strict minimizer
even though the transition lists several charges at the exact crossing.

Samples distinguish:

- `minimizing_Q` and `magnetization`: exact minimizers at the normalized sample h;
- `near_minimizing_Q` and `near_magnetization`: sectors within `tie_atol` of the
  minimum **total field-dependent energy**, including every strict minimizer;
- `status`: `tie` for multiple strict minima, `near_tie` for additional nearby
  sectors, or `unique` otherwise.

`tie_atol` affects only the nearby lists and sample status. It is not a DMRG
energy-error estimate, a transition-field tolerance, or a rule for selecting one
magnetization at degeneracy. Intervals and transitions do not depend on it.

The envelope uses O(S²) pair comparisons for S supplied sectors; sample
evaluation uses O(S) comparisons per field. Enumerating missing charges costs
O(N) time and up to O(N) output storage. Exact rational arithmetic has additional
bit-arithmetic cost. This implementation is intended for finite sector tables,
not arbitrarily large N with only a handful of energies.

## Precision labels and provenance

Optional `precision_status` is a dictionary with exactly the supplied Q keys.
Allowed strings are `not_assessed`, `not_run`, `not_evaluated`, `incomplete`,
`unmet`, and `passed`. Without the argument, all sectors are `not_assessed`.
Labels are copied from the caller and do not filter competitors or certify
energies; even a `passed` label leaves `numerical_quality="not_assessed"`.

Bare numbers cannot validate a shared model, geometry, source checkout,
environment, or numerical accuracy. The input table and output report are not
sealed run records. This unit does not import older checkpoints or bypass their
strict source/environment checks. Aggregation of compatible
[`run_static` records](static_workflow.md) and a shared Q-solve driver are future
work. Keep finite-system steps, edge rearrangements, bulk plateaus, and neutral
gaps as separate claims.

## TOML/CSV example

Run from the repository root with installed dependencies and an existing output
parent directory. The destination itself must be new:

```sh
julia --project=research --startup-file=no --threads=1 examples/magnetization_curve.jl /tmp/kagome-magnetization-demo examples/configs/magnetization_curve_demo.toml
```

The committed demo uses artificial N=4 energies `(Q,E0)=(-4,4),(-2,1),(0,0),
(2,1),(4,4)`. They are an envelope fixture, **not kagome energies**. Their
transitions occur at `h=-3,-1,1,3`. No solver or research scan runs.

The example input has an allowlisted schema:

| Key | Requirement |
| --- | --- |
| `schema_version` | Integer `1` |
| `format` | `"KagomeDMRG.magnetization_energy_table"` |
| `scope` | `"synthetic_fixture"` or `"provided_energies"`; caller-supplied description, not validation |
| `N` | Positive integer |
| `fields` | Array of sample fields; may be empty |
| `tie_atol` | Optional nonnegative absolute total-energy tolerance; default `1e-10` |
| `[[sectors]]` | Repeated table with `Q`, `energy`, and optional `precision_status` |

Unknown keys and malformed inputs are rejected before the output directory is
created. If any sector has a precision label, the example supplies the complete
status map to the API and fills omitted labels with `not_assessed`; otherwise
`precision_status_source="not_provided"` is retained.

The new directory contains `report.toml` and `samples.csv`. The TOML wrapper has
format `KagomeDMRG.magnetization_example`, schema 1, the caller's `input_scope`,
and the full public report under `curve`. CSV rows match the requested fields;
list-valued cells are quoted with semicolon-separated entries, preserving all
ties. For transitions and coverage, read the TOML report. Empty fields produce
a CSV header and no data rows.

Existing output paths, including symlinks, are rejected. Files are rendered
before directory creation, but publication of the two files is not a single
transaction: a disk error can leave a partial new directory. Retry into a new
location, preserving the earlier output for inspection.

## Validation (2026-09-23)

The focused `magnetization` group passed on Julia 1.13.0 with one Julia/BLAS
thread in the fixed research environment:

```sh
julia --project=research --startup-file=no --threads=1 -e 'using Pkg; Pkg.test("KagomeDMRG"; allow_reresolve=false, test_args=["magnetization"])'
```

The final `Pkg.test` call, including the TOML/CSV example checks, completed
successfully in 10.216 seconds excluding outer Julia startup (Pkg status output
was suppressed for that timed run). Before adding the example checks, the
direct group run passed 138 assertions in 4.7 seconds of Test time. These are
focused checks, not a new timing for the complete package suite.

Hand-calculated tables check the ordinary staircase, skipped and point-only
sectors, multiway coexistence, narrow positive intervals, incomplete coverage,
precision-label preservation, and absolute near-tie tolerance. Regression cases
cover cancellation under large common offsets, rounded crossings, collapsed
Float64 bounds, overflow, invalid input, and report serialization. The example
checks preserve every CSV tie candidate, reject invalid inputs before writing,
and retain existing output bytes when reuse is refused.

The independent N9 reference uses Cartesian-distance bonds and the spin basis,
without the production MPO builder. Ground energies from all ten Q sectors
(maximum dimension 126) are aggregated and compared with full 512-dimensional
`H0-h*sum(Sz)` diagonalization at `h=-0.2, 0.2, 4.0`. Minimum energies and
magnetizations agree within the existing small-system tolerances. No DMRG is
used in these tests.

The initial run stopped with signal 4 on the invalid `[0=>true]` input under
Julia 1.13. Validating energy before the charge calculation removed the observed
trap while preserving `ArgumentError` rejection. The same malformed cases remain
in the passing test group; charge conventions and numerical tolerances were not
changed. Independent review checked the equations, interval topology, arithmetic
limits, and metadata interpretation.

The complete package suite and Julia 1.12 were not rerun for this pure analysis
unit. No historical research fixture, additional research DMRG, N54/N72 run, or
phase-identification calculation was performed. Saved-run aggregation and
Q-by-Q optimization remain separate future work.
