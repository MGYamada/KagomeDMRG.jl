# Agent Instructions for KagomeDMRG.jl

These instructions apply throughout this repository, including delegated work.

## Research Scope and Starting Points

- Study the spin-1/2 kagome antiferromagnet at `M/Msat = 1/9` using U(1)
  cylinder DMRG and spin-flux insertion. The default target is the
  nearest-neighbor isotropic Heisenberg model.
- Read [ROADMAP.md](ROADMAP.md) and the
  [flux-insertion design](docs/flux_insertion_design.md) before changing
  physics conventions, numerical algorithms, or research protocols.
  Keep these documents consistent with implemented behavior.
- Treat the SU(3)_1 and Hall-active/Hall-inactive D(Z_3) SETs as candidate
  theories, not established phases of the microscopic model. Retain VBC,
  magnetic order, gapless states, and other competing explanations.
- Distinguish proposed work, implemented features, numerical evidence,
  converged observations, and phase-identification claims. Check the actual
  repository before assuming a planned file, API, or test already exists.
- Verify novelty claims against primary literature. Separate published
  results, preprints, research reports, and our own inferences; match model,
  geometry, charge, and flux conventions before comparing results.

## Current Research Priority

- Follow [ROADMAP.md](ROADMAP.md) for current priorities. Prioritize competing
  central bond/spin structures at nearest-neighbor 1/9 magnetization, separating
  preparation, edge, and finite-width effects. The initial N54, Q6 comparison
  of hourglass/windmill bond preparation and an unpinned complex reference is
  implemented and run, but all branches remain unconverged. The matched-parent
  chi128/256 comparison from the saved eight-sweep random and windmill states
  is also implemented and run: all four children fail all five precision
  conditions, while aligned central bond/spin/correlation differences shrink.
  Next add two fixed-chi256 sweeps to each saved ten-sweep chi256 child, as one
  bounded two-case comparison, to isolate relaxation after increasing chi.
  This next comparison is not yet implemented. Retain both branches, originally
  chosen for central bond contrast after circumferential alignment, not trial
  energy ranking. See docs/research/p4_vbc54_matched_parent_comparison.md.
  The hourglass preparation selects only part of the published
  bond hierarchy; neither it nor existing period9/27 seeds are published VBC
  wavefunctions. See the linked report in ROADMAP for the actual protocol.
- Add sectors, sweeps, bond dimensions, or geometries when they can resolve a
  stated physical uncertainty. Do not default to further N27 refinement or
  repeat passed benchmarks. Allocate resources to a physical comparison as a
  whole; previous 600/900-second run caps are not permanent research limits.
- Preserve numerical integrity and existing accuracy criteria. Unconverged
  trials can inform exploration, but do not establish converged energy rankings,
  boundaries, or phases. Use the published `h/J ≈ 0.35–0.42` only as a reference.
  At fixed Q a uniform field shifts energy by `-h*Q/2`; do not repeat DMRG over h.
- Advance the known-CSL measurement control alongside this static work.
  Its nonzero-pump calibration is required before interpreting the 1/9 pump,
  but does not block static nearest-neighbor calculations or bounded numerical
  continuation checks. Keep the control's `J2=J3=0.5,Q=0` results distinct.
- Include geometries compatible with both nine-site and 27-site competing
  order patterns. Treat edge charge rearrangements, finite-cluster field
  intervals, bulk plateaus, and neutral gaps as separate findings.
- Use [the static research protocol](docs/research/p4_static_plateau_strategy.md)
  for conventions and accuracy criteria, and the linked research reports for
  historical evidence. Do not require every combination of size, charge, seed,
  and bond dimension. Keep current decisions in ROADMAP and run history in reports.

## Active Use of Sub-agents

- **Actively deploy sub-agents for nontrivial research and development work.**
  This is a standing request to delegate: identify independent, bounded
  subtasks and launch them early without waiting for another user request.
  The lead agent should continue useful work on integration or the critical
  path while delegated work runs.
- Good parallel assignments include theory/literature checks, lattice and
  gauge audits, independent exact-diagonalization references, separate
  backend components, observables, and numerical/code review. Match the
  number of agents to useful work and available concurrency; keep trivial
  edits and inherently sequential work local.
- Give each agent a concrete objective, relevant conventions, owned files or
  a read-only scope, acceptance criteria, a computation budget, and expected
  deliverables. Ask for findings, changes, checks actually run, and unresolved
  issues rather than unfiltered logs.
- Assign one writer per file at a time. Coordinate shared interfaces,
  dependency files, fixtures, and checkpoint formats before parallel edits.
  Preserve user changes and other agents' unfinished work.
- For changes to the Hamiltonian, charge conventions, truncation, branch
  tracking, or pump measurement, use a separate agent for independent review
  when available. Reviews must check equations and numerical evidence as
  well as code structure.
- The lead agent owns resource allocation and final integration: coordinate
  CPU/BLAS threads, memory, GPU use, and separate output locations. Avoid
  duplicate expensive runs and oversubscription.
- **Keep the theta points of one flux-continuation trajectory sequential.**
  Each accepted state seeds the next point. Independent sizes, initial
  states, chirality branches, and validation runs can be delegated in parallel.
- Inspect delegated results, resolve disagreements using evidence, and run
  the necessary checks on the integrated change. An agent's completion report
  alone does not establish correctness. If delegation is unavailable,
  continue locally and state any resulting verification limitation.

## Physical Conventions

- Use `hbar = 1`. A physical `S+` raises `Sz` by 1; the integer tensor charge
  is `q = 2Sz`. For the 1/9 target, require `N % 9 == 0` and use
  `M = N/18`, `Q = N/9`, `Nup = 5N/9`, and `Ndown = 4N/9`.
  ITensor's integer `Sz` QN is `2Sz`, not physical `Sz`.
- Use open boundaries along the cylinder axis and periodic boundaries around
  its circumference. Specify the wrap vector, terminations, site ordering,
  and geometric cuts; do not infer a literature YC label from `Ly` alone.
- Store oriented bonds with signed winding `wy`. In seam gauge, the
  coefficient of `S+_i S-_j` is `(Jxy/2) * cis(wy * theta)` and the reverse
  term has its complex-conjugate coefficient. Do not twist `Sz_i Sz_j`.
  Count each physical bond once and reverse winding with bond orientation.
- Retain complex states and Hermitian arithmetic, including at zero flux
  when exploring complex branches. Distinguish transpose from adjoint.
  Explicit symmetry-breaking seeds or added couplings change the model:
  label them and check their removal before claiming a nearest-neighbor result.
- Measure spin transfer across a cut from the sum of
  `Sz(theta) - Sz(0)` in the right region, using the measured zero-flux
  profile. Preserve unwrapped theta and the full cumulative transfer.
  Do not round the transfer modulo 1 or identify `dE/dtheta` with axial pumping.
- The dimensionless response is `sigma_xy^Sz = t' * inv(K) * t` in the
  conventions of the design document. Expected values such as `2/3` and `0`
  are hypotheses, never branch-selection or convergence criteria.
- A `2/3` pump alone does not distinguish SU(3)_1 from Hall-active D(Z_3),
  and a zero pump alone does not establish a trivial phase. For the documented
  candidates, the six-pi cycle restores the bulk topological sector, not
  necessarily the full finite-cylinder wavefunction or its edge densities.
  A three-cycle flux orbit does not establish threefold total degeneracy.

## Numerical Validation and Research Claims

- Validate bond counting, Hermiticity, U(1) conservation, winding signs,
  zero flux on contractible loops, circumference flux, seam-gauge two-pi
  periodicity, and gauge equivalence before interpreting a pump.
- Compare small systems at zero and generic nonzero theta against an
  independent spin-basis exact-diagonalization implementation. Keep reference
  construction independent enough to expose shared bond/phase errors; do not
  merely call the production MPO builder from the reference test.
- Validate the measurement pipeline with a known CSL pump and a controlled
  zero-response case before interpreting the 1/9 spin pump or Hall response. Do not demand
  fractional pumping from a small system with a unique periodic ground state.
- Warm starts and low energy do not guarantee branch continuity. Combine
  overlaps, bulk density, entanglement/Schmidt charge, residuals or variance,
  and measured truncation errors. On discontinuity, restore the last accepted
  checkpoint and refine the step; mark unresolved trajectories explicitly.
- Check that left and right transfers cancel, that several bulk cuts agree,
  and that real-space transfer matches calibrated Schmidt-charge changes.
  Examine bond-dimension, flux-step, length, width, and initial-state dependence
  before claiming quantization. A configured cutoff is not an error estimate.
- Include geometries compatible with competing enlarged unit cells. Separate
  edge rearrangements and domain-wall motion from bulk transport. A magnetization
  plateau does not prove a neutral bulk gap; identify edge-localized excitations.
- Use additional sector, entanglement, or modular diagnostics for intrinsic
  topological order. Do not infer chiral central charge from one entanglement
  spectrum or a nonzero chirality expectation alone.

## Implementation and Verification Workflow

- Establish an ITensors.jl + ITensorMPS.jl reference backend before replacing
  numerical kernels. Share explicit lattice, model, observable, and checkpoint
  contracts across backends. Compare custom-backend results with both ED and
  the reference at matched accuracy.
- Use `../SUNDMRG.jl` as a design reference for sweeps, workspaces, storage,
  and parallel execution. Its SU(N) irreps, multiplet weights, coefficients,
  and real-valued kernels are not a drop-in U(1) flux implementation.
- Preserve charge conservation and global truncation weights across sectors.
  Allow the bond charge sectors needed for transport to grow; do not freeze
  their set during continuation. Refresh Hamiltonian-dependent environments
  after a flux update, or maintain correctly updated separate components.
- Profile before adding GPU/MPI or rewriting tensor storage. Compare elapsed
  time and memory at the same physical accuracy, including transfer and I/O
  overhead. Verify complex block-sparse behavior on each backend used.
- Run focused, independent numerical checks appropriate to the change and
  the existing test suite when applicable. Once the Julia project and test
  entry point exist, the standard package check is:

  ```sh
  julia --project=. --startup-file=no -e 'using Pkg; Pkg.test()'
  ```

- Do not launch large DMRG scans or set up an accelerator environment merely
  to validate a documentation-only change. Report checks actually performed,
  failures, and remaining limitations without claiming unrun tests passed.

### Practical Test-Time Budget

- **Keep routine testing fast enough for daily research work.** With dependencies
  installed, aim for focused checks within about 60 seconds. Keep the complete
  package suite within about 3 minutes when practical; 5 minutes is a trigger to
  investigate and reduce cost, not a reason to silently omit failures. Record
  Julia startup/precompilation separately when it dominates elapsed time.
- Select affected groups with `Pkg.test(; test_args=["dmrg", "checkpoint"])`
  or the corresponding groups listed in `test/README.md`. Run the complete suite
  at most once for an integrated change that warrants it. Repeat only checks
  affected by a subsequent fix; do not repeat passed numerical work after edits
  to prose, labels, or reports.
- Prefer a small set of cases that expose distinct failures. Delete redundant
  cases and Cartesian products rather than retaining them in an exhaustive
  mode. Keep known-bug regressions and independent numerical references, and
  never relax physical tolerances merely to meet a time budget.
- If validation exceeds the budget, identify the expensive case or compilation
  step, reuse bounded fixtures, and narrow further runs. Record checks omitted
  and the reason; do not add more test infrastructure or tests just to increase
  assertion counts. Count reduction alone does not establish a speedup.
- Keep research scans, convergence studies, and repeated fresh-process checks
  outside the routine edit/test loop. Before running them, set an explicit
  question, cases, CPU/thread allocation, and computation bound. Save their
  outcomes as research evidence rather than turning every point into a package
  regression. Necessary checkpoint/provenance regressions remain in the suite.

## Reproducibility and Saved Results

- Save explicitly selected metadata: model/couplings, geometry and ordering,
  total charge, gauge and signs, seed, accepted flux path, solver settings,
  code revision, dependency versions, and numerical diagnostics.
- Keep accepted checkpoints separate from trial states. Save atomically and
  verify geometry, site indices, charge, and configuration before resuming.
  Use distinct output locations for independent trajectories.
- Retain zero, nonquantized, failed, and competing-state outcomes with their
  status and uncertainty. Keep exploratory results separate from converged
  claims; update the roadmap when milestones are actually achieved.
- Use an allowlist for saved metadata. Do not dump all environment variables,
  host settings, or connection configuration for reproducibility.

## Local Network Configuration

- Do not record the user's local network configuration or connection details in
  repository files, documentation, source comments, test fixtures, generated
  artifacts, saved diagnostic logs, commit messages, issue/PR descriptions,
  or persistent agent notes.
- This includes private IP addresses, subnet/gateway/DNS settings, local
  hostnames, SSIDs, router settings, port forwarding and firewall rules,
  and SSH endpoints, usernames, key paths, and credentials.
- Use connection details only as needed for the current authorized task.
  Use placeholders such as `<GPU_HOST>` and `<SSH_USER>` in saved examples,
  launch commands, and instructions.
- Redact network and connection details before saving diagnostic output or
  publishing validation reports.
- GPU/MPI validation records may include GPU models/counts, OS and software
  versions, numerical results, and commands with placeholders. Omit local
  network setup and remote-access configuration.
