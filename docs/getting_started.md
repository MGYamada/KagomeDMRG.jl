# Getting started

KagomeDMRG.jl is a research package for spin-1/2 kagome cylinders, with the
nearest-neighbor antiferromagnetic Heisenberg model at `M/Msat = 1/9` as its
primary target. It currently provides lattice and gauge construction, a
complex U(1) two-site DMRG reference solver using ITensors.jl and ITensorMPS.jl,
and basic observables. Completed validation and its dependency baselines are
recorded below; larger research systems have not been validated.

The current priority is to compare competing central bond/spin structures on
the same longer cylinder, separating preparation, boundary, and finite-width
effects. The planned N54 comparison uses hourglass/windmill-inspired preparation
and an unpinned complex reference; that protocol is not yet implemented.
Additional sectors and higher accuracy are chosen to resolve specific physical
uncertainties. Known-CSL pump calibration proceeds alongside static research
and remains required before interpreting a 1/9 quantized pump. See the
[roadmap](../ROADMAP.md) for priorities and the
[static protocol](research/p4_static_plateau_strategy.md) for conventions and criteria.

## Run a small system

Run these commands from the repository root with Juliaup installed. The
directory override selects Julia 1.13 for this checkout; the explicit
`+1.13` commands also work without an override.

```sh
juliaup add 1.13
juliaup override set 1.13
julia +1.13 --project=research --startup-file=no -e 'using Pkg; Pkg.instantiate(); Pkg.precompile()'
julia +1.13 --project=research --startup-file=no -e 'using Pkg; Pkg.test("KagomeDMRG")'
julia +1.13 --project=research --startup-file=no --threads=1
```

`research/Project.toml` is the dedicated research environment.
`research/Manifest-v1.13.toml` records its Julia 1.13 dependencies. Julia
automatically selects the manifest matching its major and minor version,
as described in the [Pkg documentation](https://julialang.github.io/Pkg.jl/v1/toml-files/#Different-Manifests-for-Different-Julia-versions).
`research/Manifest.toml` records its Julia 1.12.7 dependencies.
Both manifests select this checkout and the local NDTensors correction by relative
paths. Use `--project=research` for the research examples. Historical validation
reports record their original source and manifest hashes; reproducing those
older results requires the corresponding historical checkout.
The compatibility ranges in `Project.toml` already allow Julia 1.13.
To select Julia 1.12.7, run `juliaup add 1.12.7` and use
`julia +1.12.7` in the Julia commands instead.

In the research environment, `Pkg.test("KagomeDMRG")` runs the regression suite.
Use `Pkg.test("KagomeDMRG"; test_args=["truncation", "checkpoint"])` to select groups.
The library root remains a separate development workflow:
`julia --project=. --startup-file=no -e 'using Pkg; Pkg.instantiate(); Pkg.test()'`.
It resolves the root project's compatibility ranges and may create an untracked
root `Manifest*.toml`; use the manifests in `research/` for the committed research
baseline. Existing checkpoints retain their original source/environment identity
and require the corresponding historical checkout for restart.
See the [test guide](../test/README.md) for the available groups and coverage.

At the Julia prompt:

```julia
using KagomeDMRG, LinearAlgebra, ITensors
BLAS.set_num_threads(1)
ITensors.disable_threaded_blocksparse()

lattice = kagome_cylinder(1, 3)
sector = target_sector(nsites(lattice))
result = run_dmrg(lattice, 0.37; seed=11)

@show sector                       # (N=9, Q=1, M=0.5, Nup=5, Ndown=4)
@show result.energy                # expectation value of the final Hamiltonian
@show result.sz sum(result.sz)     # physical local and total Sz
@show result.max_truncation_errors # measured maximum for each sweep
@show result.variance              # raw <H^2> - <H>^2
```

Flux `theta` is in radians and remains unwrapped. `run_dmrg` optimizes one
flux point; it does not certify convergence or branch continuity. Its default
settings are for small validation systems. The variance may be slightly
negative from floating-point cancellation. Measured truncation errors are
separate from the configured `cutoff` and do not by themselves bound the
error of a physical observable.

The [validation example](../examples/validate_small_system.jl) repeats the
small-system ED comparisons and saves selected metadata and numerical metrics:

```sh
julia +1.13 --project=research --startup-file=no examples/validate_small_system.jl
```

An optional positional argument selects a new output TOML file. Use a distinct
output location for each run and retain failed or nonconverged outcomes.

## Geometry, charge, and gauge conventions

Cell coordinates are `x=0:Lx-1` with open boundaries and `y=0:Ly-1` with
periodic boundaries. The primitive vectors are `a1=(1,0)` and
`a2=(1/2,sqrt(3)/2)`, and sublattice offsets are `A=0`, `B=a1/2`, and `C=a2/2`.
The cylinder wrap vector is `Ly*a2`. `Ly` counts primitive cells and does not
automatically specify a literature YC label. Only `Lx >= 1` and `Ly >= 3`
are supported.

Sites are ordered by increasing `x`, then `y`, then `A,B,C`, with one-based
site indices. The lattice keeps all sites in these cells and removes bonds
whose `x` endpoint lies outside the cylinder. It has `N=3Lx*Ly` sites and
`(6Lx-2)*Ly` nearest-neighbor bonds. `site_index` requires canonical coordinates
and does not wrap out-of-range inputs.

With `hbar=1`, a physical `S+` raises `Sz` by one. The integer tensor charge
is `q=2Sz`, so the target sector requires `N % 9 == 0` and uses
`Q=N/9`, `M=N/18`, `Nup=5N/9`, and `Ndown=4N/9`. The nine-site example thus
has physical total `Sz=1/2`, not `Sz=1`.

For a different fixed sector, pass an integer `Q` to `target_sector`,
`initial_mps`, or `run_dmrg`, for example
`run_dmrg(kagome_cylinder(1, 4), 0.37; Q=0, seed=11)`.
Require `-N <= Q <= N` and matching parity of `N` and `Q`; the up-spin count
is `(N+Q)/2`. Omitting `Q` still selects 1/9, even with a supplied `psi0`:
a warm start with another charge must name that charge explicitly.

The bounded [charge-sector study](research/p3_explicit_charge_validation.md)
compares twelve-site `Q=0` and nine-site neighboring sectors against independent
ED, outside the regression suite:

```sh
julia --project=research --startup-file=no --threads=1 examples/validate_charge_sectors.jl outputs/new-charge-study
```

Choose a new output directory. The study records all five points and their
accuracy, with a 300-second budget checked between points. It does not measure
an axial pump or establish a plateau.

Each physical bond is stored once as `(i,j,Jxy,Jz,wy)`, with `wy` the signed
periodic image of endpoint `j` relative to `i`. Reversing the bond reverses
`wy`. The coefficient of `S+_i S-_j` is `(Jxy/2)*cis(A_b)`, its reverse is the
complex conjugate, and `Sz_i Sz_j` is not twisted. In seam gauge,
`A_b=wy*theta`. In uniform gauge, add `chi_i-chi_j`, where
`chi_i=-theta*eta_i/Ly` and `eta_i=y_i+(sublattice_i == :C ? 1/2 : 0)`.
`bond_phase` returns the real angle `A_b`, not its exponential.

`right_region(lattice, cut)` selects every site in cells with `x >= cut`,
with `1 <= cut < Lx`. These are cuts between complete cell columns. The
`Lx=1` quickstart has no such interior cut. For larger systems,
`spin_transfer(lattice, baseline, current)` sums `current-baseline` in the
left and right regions of each cut. Supply the measured zero-flux profile as
`baseline`; transfer is kept cumulative and is never rounded modulo one.

## Available API

| Purpose | Functions |
| --- | --- |
| Geometry and sectors | `kagome_cylinder`, `nsites`, `site_index`, `target_sector`, `right_region` |
| Extended isotropic exchange | `kagome_j1j2j3_cylinder`, `bond_families` |
| Bond orientation and gauge | `reverse_bond`, `bond_phase`, `gauge_angles` |
| U(1) reference backend | `spin_sites`, `initial_mps`, `twisted_exchange_mpo`, `run_dmrg` |
| Observables | `sz_profile`, `spin_correlations`, `spin_transfer`, `bond_energies` |
| Schmidt probabilities and absolute left charge | `schmidt_diagnostics` |
| Local snapshots and completed-point restart | `save_checkpoint`, `load_checkpoint`, `resume_dmrg` |
| Diagnostic-gated adaptive continuation | `FluxPolicy`, `continue_flux` |

The exported data types are `Bond`, `KagomeSite`, `KagomeCylinder`, and `FluxPolicy`.
Use Julia help, for example `?run_dmrg`, for arguments and result fields.

`kagome_j1j2j3_cylinder(Lx,Ly; J1=1.0,J2=0.5,J3=0.5)` constructs the
separate extended model used to prepare a CSL control. J3 includes only
opposite hexagon vertices. A bond remains when both endpoints survive the
open boundary, even if the rest of its hexagon does not; zero-coupling bonds
also remain explicitly. The NN builder's defaults are unchanged.
`bond_families(lattice)` classifies actual bonds as `:J1`, `:J2`, `:J3`, or
`:other` from their geometry and winding; checkpoints save and validate these
families together with every exchange coefficient.

`bond_energies(psi,lattice,theta; gauge=:seam)` returns exchange-energy
expectations in bond order. Fields are excluded: subtract
`dot(hz,sz_profile(psi))` from their sum to obtain the full energy.
This implementation uses the full correlation matrices and has O(N²) output
storage; large-cylinder cost has not been profiled.
See the [extended-model study](research/p3_extended_model_validation.md) for
the independent small-system validation scope.

```julia
control = kagome_j1j2j3_cylinder(1, 4; J1=1.0, J2=0.5, J3=0.5)
point = run_dmrg(control, 0.37; Q=0, seed=11)
exchange = bond_energies(point.psi, control, point.theta)
```

Reproduce the bounded two-point study with a new output directory:

```sh
julia --project=research --startup-file=no --threads=1 examples/validate_extended_model.jl outputs/new-extended-study
```

`spin_correlations` returns `zz=<Sz_i Sz_j>` and `pm=<S+_i S-_j>`; transverse
correlations can be complex. `run_dmrg` returns the MPS and MPO as `psi` and
`H`, the site indices, flux, gauge, integer charge, observables, per-sweep
diagnostics, and solver settings.

To supply an initial state, create `sites=spin_sites(lattice)` once and pass
the same indices to `initial_mps`, `twisted_exchange_mpo`, and `run_dmrg`.
An optional `psi0` must use those exact indices and the target total charge;
the solver copies it before optimization. A warm start does not establish
that a physical branch was followed.

## Checkpoint restart and Schmidt diagnostics

The checkpoint API saves a completed optimization point together with its
measured zero-flux baseline. For example, with the nine-site `lattice` above:

```julia
baseline = run_dmrg(lattice, 0.0; seed=11)
point = run_dmrg(lattice, 0.37; sites=baseline.sites, psi0=baseline.psi, seed=11)
snapshot = save_checkpoint("outputs/restart-demo", point;
    baseline, theta_path=[0.0, 0.37], status=:accepted)
restored = load_checkpoint(snapshot, lattice;
    expected_theta=0.37, expected_settings=point.settings)
next_point = resume_dmrg(snapshot, lattice, 0.71;
    expected_settings=point.settings)
schmidt = schmidt_diagnostics(next_point.psi, 4)
@show schmidt.entropy schmidt.mean_left_sz
```

`save_checkpoint` defaults to `status=:trial`. Trial and accepted snapshots
occupy separate directories, and each save creates a new immutable directory.
The caller supplies the path and acceptance status; saving an accepted snapshot
does not establish convergence or branch continuity. `resume_dmrg` accepts only
accepted snapshots and starts a new sweep batch using their settings and exact
site indices, rebuilding the Hamiltonian at the requested unwrapped flux.
It does not resume inside an unfinished sweep or advance an adaptive trajectory.
`load_checkpoint`, `resume_dmrg`, and `continue_flux` inherit the saved or
starting state's validated `Q` when it is omitted. An explicit `Q` is an
expectation and must match; these functions never change the charge mid-run.
Schema 1 already stores the charge in model and state metadata, so the schema
is unchanged. Snapshots from older source versions remain incompatible under
the strict source-identity check; there is no automatic migration.

Each snapshot contains the current and zero-flux MPS, allowlisted TOML metadata,
and byte lengths and SHA-256 checksums. Loading checks the model, full geometry,
bonds, field, gauge, charge, site identity, solver settings, baseline, state
normalization, and energy. Julia, package versions, architecture, source hashes,
and the active Project plus Julia's selected dependency manifest must match.
The primitive metadata and integrity envelope are checked before the Julia
payload is deserialized. Package and vendored source hashes are captured when
modules are evaluated and remain precompile dependencies. The active environment
and runtime are captured afresh in `__init__`, even when compiled code is reused.
The result and baseline carry that immutable execution identity. Editing source,
Project, or Manifest files, or switching the active project, invalidates later
runs, saves, loads, and continuation in that process; restart Julia after such
changes. Saved environment hashes use role-prefixed filenames, without absolute
paths. Saving never relabels an older in-memory result with a new identity.
Use these snapshots only from trusted local runs. They are deliberately bound
to the same environment, not a portable archival format; see the
[Julia Serialization contract](https://docs.julialang.org/en/v1/stdlib/Serialization/).
Publishing a fully written directory uses a same-parent rename. This provides
atomic visibility, not a guarantee against power loss. Saving the baseline MPS
in every snapshot costs additional storage; production-scale I/O is unprofiled.

`schmidt_diagnostics(psi, b)` measures the MPS prefix `1:b`, preserves its input,
and performs a block-sparse decomposition without truncation. It returns
probabilities, one absolute integer `left_q=2Sz_left` per Schmidt state,
natural-log entropy, mean and variance of physical left `Sz`, and input norm².
The charge origin and sign come from tensor flux conservation, independently
of the density profile. For a geometric cylinder cut `c`, use `b=3lattice.Ly*c`.
The nine-site example has no interior geometric cut: bond 4 is an MPS cut only.

The separate-process restart and independent ED checks can be reproduced with:

```sh
julia +1.13 --project=research --startup-file=no --threads=1 examples/validate_restart.jl
```

This writes snapshots and an allowlisted numerical report to a new
`outputs/p1-restart-*` directory. The manually selected flux points test restart
semantics, not an adiabatic flux trajectory. The caller may supply a new output
directory as the only positional argument.

## Adaptive flux continuation

`continue_flux` advances through an ordered list of unwrapped target angles.
It reloads the last accepted checkpoint before every trial, rebuilds the MPO,
and halves the step after a failed continuity or convergence check. A reduced
step persists after acceptance; the driver does not raise bond dimension or
sweep count automatically. All thresholds in `FluxPolicy` must be supplied
explicitly except `consistency_tol`. They must be appropriate to the system
size, monitored region and solver accuracy; no expected pump value is used.

This executable control has a unique, flux-independent product ground state.
It changes the Hamiltonian to longitudinal fields only and tests the pipeline,
not the phase of the Heisenberg model:

```julia
using KagomeDMRG, ITensors, ITensorMPS, LinearAlgebra
BLAS.set_num_threads(1)
ITensors.disable_threaded_blocksparse()
control = kagome_cylinder(3, 3; Jxy=0, Jz=0)
sites = spin_sites(control)
labels = [fill("Up", 15); fill("Dn", 12)]
hz = [fill(1.0, 15); fill(-1.0, 12)]
start = run_dmrg(control, 0.0; sites, psi0=MPS(ComplexF64, sites, labels),
    hz, nsweeps=2, maxdim=2, cutoff=0.0)
policy = FluxPolicy(min_overlap=1-1e-10, max_density_change=1e-10,
    max_entropy_change=1e-10, max_schmidt_change=1e-10,
    max_variance=1e-9, max_truncation_error=1e-10,
    max_sweep_energy_change=1e-9, max_cut_spread=1e-9)
trajectory = continue_flux(control, [2pi, 4pi, 6pi, 0.0]; start,
    output_root="outputs/zero-control", policy, initial_step=2pi, min_step=pi/32,
    hz, cuts=[1,2], diagnostic_bonds=[9,18], bulk_sites=10:18)
@show trajectory.status trajectory.theta_path trajectory.output_path
```

The first version requires seam gauge, zero DMRG noise, measured variance at
every point, and at least two sweeps. `min_overlap` refers to overlap amplitude,
not squared fidelity. Other successive-point checks compare monitored-site
density, entropy and mean left Schmidt spin. Quality checks use the last
sweep's measured truncation error, variance, the difference of the last two
local sweep energies, and their difference from the measured final energy.
Negative variance must lie within `100eps(Float64)*max(1,E^2)` and the explicit
variance limit. Every sweep's raw diagnostics remain in the record.

Transfer is always measured against the original zero-flux state. The driver
checks left/right cancellation and agreement with calibrated Schmidt changes,
and records the full transfer and its spread across selected cuts. A cut-spread
failure can be physical bulk rearrangement and does not imply a software error.
Every geometric cut is added to the diagnostic Schmidt bonds. `bulk_sites`
defaults to all sites; choose and record actual interior sites for a larger
cylinder. An auxiliary Schmidt bond on the nine-site cylinder does not create
an axial geometric cut.

A run creates a unique directory with an atomically updated `trajectory.toml`
journal, explicit policy and selections, separate trial and accepted snapshots,
and rejection reasons. `status == :completed` means the finite diagnostic
limits passed at the targets. It is not a certificate of a physical branch,
quantization or phase identity. `:unresolved` records minimum-step or trial-budget
exhaustion, invalid diagnostics or other solver errors. Valid rejected states
remain under `trial/`; invalid states cannot be promoted. Interrupts, resource
exhaustion or an I/O failure may leave a journal marked `running`, while already
published checkpoints remain available.

To restart, pass an accepted checkpoint path as `start`, together with the
same lattice, `hz` and gauge and explicit new targets, policy and step settings.
The solver settings, original baseline and accepted path are preserved, and a
new journal is created. This is a restart at a completed point, not restoration
of an interrupted driver instruction. On an initial quality failure,
`last_checkpoint === nothing` and no point is newly accepted.

The bounded reproduction includes a restart, forward/reverse zero control,
independent nine-site ED, and a degenerate initial state that stays unresolved:

```sh
julia +1.13 --project=research --startup-file=no --threads=1 examples/validate_continuation.jl
```

See the [P2 validation report](research/p2_continuation_validation.md). The known
CSL positive control and convergence studies on longer interacting cylinders
remain pending.

## Validation and present limits

For a bounded interacting-cylinder study with an axial cut, run:

```sh
julia +1.13 --project=research --startup-file=no --threads=1 examples/validate_interacting18.jl
```

This fixes the nearest-neighbor Heisenberg model in the `N=18,Q=2` sector and
compares a lower-rank initial state with a state retaining the full fixed-charge
rank. It compares steps `0.185` and `0.0925`, seeds 11 and 29, forward/return
paths through `theta=0.37`, and negative flux. Every saved point is checked
against an independent sparse BlockLanczos reference; low-accuracy and rejected
outcomes remain in the report. The partial ED spectrum is calibrated against
dense nine-site results, including degenerate eigenspaces. Its convergence count
and full-matrix residuals are checked explicitly.

There is only one geometric cut at bond 9 and no interior bulk column in this
two-column system. The study does not establish a magnetization plateau, bulk
gap, quantized pump, or microscopic phase. All eigensolver settings and finite
diagnostic limits are saved with the selected source hashes.

The [18-site report](research/p2_interacting18_validation.md) records four
completed chi=512 paths and 24 saved-point comparisons, including reused initial
states. Every point passed the independent ED limits; halving the preset step or
changing the seed changed the measured transfer by at most `4.21e-13`.
Chi=128 failed the initial accuracy gates. A separate chi=256 run stopped at the
truncation-spectrum guard before producing a final state; this failure is also
preserved. No chi=512 trial was rejected, so these paths did not exercise adaptive
step refinement. The literature audit and implementation requirements for the
known CSL control are in the [P3 design](research/p3_csl_control_design.md).

The following numerical results predate the research-environment relocation.
Their original root manifest paths, hashes, and commands remain historical records;
moving the current environment does not revalidate or relabel those results.

On 2026-09-19, the standard `Pkg.test()` command passed 2,205 assertions with
Julia 1.12.7, ITensors 0.9.31, ITensorMPS 0.4.1, and KrylovKit 0.10.4.
This is a tested baseline, not validation of every version allowed by the
compatibility ranges. For `theta=0` with seed 11 and `theta=0.37` with seeds
11 and 29, the maximum energy error per site was `4.94e-16`, the maximum
local Sz error was `5.25e-14`, and the maximum ED residual norm was `2.31e-13`.
See the [validation report](research/p0_reference_validation.md) for the scope
of these results.

On the same date, Julia 1.13.0 also passed all 2,205 assertions using
`Manifest-v1.13.toml`, with ITensors 0.9.31, ITensorMPS 0.4.1, and KrylovKit
0.10.4. Dependency resolution, instantiation, and precompilation completed
successfully. The numerical error values above belong to the earlier
Julia 1.12.7 run. The small-system validation example also passed all three
independent points on Julia 1.13.0 and saved matching pre-run and post-run
SHA-256 hashes for `Manifest-v1.13.toml`.

After adding checkpoint and Schmidt diagnostics, Julia 1.13.0 passed all
3,782 assertions: the original 2,205, plus 1,477 Schmidt checks and 100
checkpoint checks. The Julia 1.12.7 manifest was re-resolved only to add
direct standard-library dependencies; the expanded suite has not been run
on Julia 1.12.7. See the [restart and Schmidt report](research/p1_restart_schmidt_validation.md).
The separate-process restart reproduced the direct continuation energy and
local Sz profile exactly in this run; the overlap magnitude error was
`4.44e-16`. The maximum independent ED residual over the recorded small-system
checks was below `1.06e-13`. These are completed-point restart checks, not
evidence of an adiabatic branch or quantized pumping.

After adding adaptive continuation and the truncation rank guard, the Julia
1.13.0 integrated `Pkg.test()` run passed all **4,283 assertions**: the previous
3,782 plus 417 continuation checks and 84 truncation checks. The standalone P2
reproduction passed with matching source hashes, including the expected
`unresolved` outcome for a degenerate initial state. See the
[P2 report](research/p2_continuation_validation.md) for numerical values and scope.

Adding the sparse-reference tests brought the Julia 1.13.0 integrated suite to
**4,348 passing assertions** (4,283 previous checks plus 65 independent
eigensolver, degeneracy and nonconvergence checks). The 18-site study is an
explicit research example, separate from the routinely bounded package tests.
That is a historical full-suite result. For the local truncation correction,
3,783 independent calibration assertions and 346 backend regressions passed.
After the final provenance representation change, 109 atomic-checkpoint and
20 focused source-drift/restart assertions passed. The expanded full-suite run
was deliberately interrupted during compilation; it is not reported as a full
pass. See the [validation scope of that change](research/p1_qn_truncation_calibration.md#テスト実行の範囲).
Redundant cases and combinations were subsequently removed. Current commands,
coverage, and run results are in the
[test guide](../test/README.md).

Tests compare geometry and signed periodic images with a Cartesian distance
oracle and compare the MPO with an independent spin-basis Hamiltonian at
zero and generic nonzero flux. DMRG tests compare energies, residuals,
densities, and correlations against independent exact diagonalization (ED).
For the nine-site cylinder in the target sector, zero flux has a twofold
ground space. The check projects the DMRG state onto the full ED ground space
and compares observables in that matched state; comparing with an arbitrary
single degenerate eigenvector would be inappropriate. This small-system
degeneracy is not a claim about topological degeneracy.

Completed-point checkpoint/restart, Schmidt-charge diagnostics, and adaptive
flux acceptance/refinement/rollback are implemented and tested on bounded
controls. Reproduction of a known chiral spin-liquid pump remains unverified.
A longitudinal-field product-state control tests zero response and its
forward/reverse continuation for an explicitly altered Hamiltonian.
It does not identify the phase of the nearest-neighbor model. No plateau,
quantized pump, or microscopic phase identification is claimed. SU(3)₁ and
Hall-active or Hall-inactive D(Z₃) remain candidate theories alongside
competing ordered or gapless explanations.

An analytic two-spin control checks a known nonzero truncation error and the
difference between the local eigensolver energy and the final MPS energy.
The previously tested upstream NDTensors 0.4.31 can discard every state when `maxdim` splits
equal or nearly equal Schmidt weights across sectors. It can also retain a
nonzero normalized state while underreporting its discarded probability:
weights `[0.6,0.2,0.2]` and `maxdim=2` retain only rank one, reporting 0.2 loss
instead of 0.4. The observer rejects zero/nonfinite states and inconsistent
reported-spectrum versus retained-link dimensions. The pinned local NDTensors
`0.4.31+1` now uses the same globally selected entries for the actual block ranks
and reported spectrum. `maxdim` is a hard cap; exact ties use block-coordinate,
then within-block-index order and may split a degenerate subspace. Positive
`min_blockdim` constraints consume that same global cap and incompatible
constraints throw an error. General non-Hermitian eigen retains the upstream
path. See the [vendor contract](../vendor/README.md).
Independent complex two/four-spin dense references check SVD, Hermitian eigen,
both factorization directions, norm scaling and cutoff conventions. With
noise, perturbed-density trace loss and original-wavefunction loss are checked
separately. Correct local truncation loss does not bound physical-observable
errors or establish convergence of a finite-bond-dimension optimization.
The nine-site reference runs use enough bond dimension to avoid truncation.
See the [recorded boundary case](research/data/qn_truncation_boundary.toml).
The [correction and finite-chi study](research/p1_qn_truncation_calibration.md)
records the captured 18-site failure and its replay. Both research manifests select the
local backend; loading an unrelated registry NDTensors is rejected. When using
KagomeDMRG from another project, that project must select the same local source.

Current discussion and research plans may remain in Japanese. Public API
names, docstrings, and machine-readable result keys use English. This guide
is an initial English entry point; full English research documentation and
release checks remain future work. See the [roadmap](../ROADMAP.md) and
[physical and numerical design](flux_insertion_design.md) for the research
protocol.
