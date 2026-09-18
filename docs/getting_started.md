# Getting started

KagomeDMRG.jl is a research package for spin-1/2 kagome cylinders, with the
nearest-neighbor antiferromagnetic Heisenberg model at `M/Msat = 1/9` as its
primary target. It currently provides lattice and gauge construction, a
complex U(1) two-site DMRG reference solver using ITensors.jl and ITensorMPS.jl,
and basic observables. The package test suite passes on the dependency
baseline recorded below; larger research systems have not been validated.

## Run a small system

Run these commands from the repository root. The Julia compatibility range
is specified in `Project.toml`; `Manifest.toml` records the resolved dependency
baseline.

```sh
julia --project=. --startup-file=no -e 'using Pkg; Pkg.instantiate()'
julia --project=. --startup-file=no -e 'using Pkg; Pkg.test()'
julia --project=. --startup-file=no --threads=1
```

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
julia --project=. --startup-file=no examples/validate_small_system.jl
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
| Bond orientation and gauge | `reverse_bond`, `bond_phase`, `gauge_angles` |
| U(1) reference backend | `spin_sites`, `initial_mps`, `twisted_exchange_mpo`, `run_dmrg` |
| Observables | `sz_profile`, `spin_correlations`, `spin_transfer` |

The exported data types are `Bond`, `KagomeSite`, and `KagomeCylinder`.
Use Julia help, for example `?run_dmrg`, for arguments and result fields.
`spin_correlations` returns `zz=<Sz_i Sz_j>` and `pm=<S+_i S-_j>`; transverse
correlations can be complex. `run_dmrg` returns the MPS and MPO as `psi` and
`H`, the site indices, flux, gauge, integer charge, observables, per-sweep
diagnostics, and solver settings.

To supply an initial state, create `sites=spin_sites(lattice)` once and pass
the same indices to `initial_mps`, `twisted_exchange_mpo`, and `run_dmrg`.
An optional `psi0` must use those exact indices and the target total charge;
the solver copies it before optimization. A warm start does not establish
that a physical branch was followed.

## Validation and present limits

On 2026-09-19, the standard `Pkg.test()` command passed 2,205 assertions with
Julia 1.12.7, ITensors 0.9.31, ITensorMPS 0.4.1, and KrylovKit 0.10.4.
This is a tested baseline, not validation of every version allowed by the
compatibility ranges. For `theta=0` with seed 11 and `theta=0.37` with seeds
11 and 29, the maximum energy error per site was `4.94e-16`, the maximum
local Sz error was `5.25e-14`, and the maximum ED residual norm was `2.31e-13`.
See the [validation report](research/p0_reference_validation.md) for the scope
of these results.

Tests compare geometry and signed periodic images with a Cartesian distance
oracle and compare the MPO with an independent spin-basis Hamiltonian at
zero and generic nonzero flux. DMRG tests compare energies, residuals,
densities, and correlations against independent exact diagonalization (ED).
For the nine-site cylinder in the target sector, zero flux has a twofold
ground space. The check projects the DMRG state onto the full ED ground space
and compares observables in that matched state; comparing with an arbitrary
single degenerate eigenvector would be inappropriate. This small-system
degeneracy is not a claim about topological degeneracy.

P1 checkpoint/resume support is still pending. Adaptive flux continuation,
rollback, Schmidt-charge calibration, and reproduction of a known chiral
spin-liquid pump are not implemented. A longitudinal-field product-state
control tests zero-response readout for an explicitly altered Hamiltonian.
It does not identify the phase of the nearest-neighbor model. No plateau,
quantized pump, or microscopic phase identification is claimed. SU(3)₁ and
Hall-active or Hall-inactive D(Z₃) remain candidate theories alongside
competing ordered or gapless explanations.

An analytic two-spin control checks a known nonzero truncation error and the
difference between the local eigensolver energy and the final MPS energy.
The tested upstream QN backend can discard every state when `maxdim` splits
exactly degenerate Schmidt weights across sectors. The observer rejects zero
or nonfinite states immediately. This is a failure guard, not an upstream
kernel fix; general degenerate-boundary truncation errors remain unvalidated.
The nine-site reference runs use enough bond dimension to avoid truncation.
See the [recorded boundary case](research/data/qn_truncation_boundary.toml).

Current discussion and research plans may remain in Japanese. Public API
names, docstrings, and machine-readable result keys use English. This guide
is an initial English entry point; full English research documentation and
release checks remain future work. See the [roadmap](../ROADMAP.md) and
[physical and numerical design](flux_insertion_design.md) for the research
protocol.
