# Explicit density/bond-biased fixed-charge trial states. These motifs are
# not wavefunctions or bond patterns taken from a published VBC proposal.
using LinearAlgebra
using ITensors
using ITensorMPS
using KagomeDMRG

function _order_seed_cell(kind::Symbol, x::Integer, y::Integer, origin)
    u, v = x - origin[1], y - origin[2]
    minus = mod(u - v, 3)
    pair_index = kind === :period9 ? minus : mod(u + v, 3)
    pair = ((1, 2), (2, 3), (1, 3))[pair_index + 1]
    return (; pair, spectator=6-sum(pair), charge=minus == 2 ? -1 : 1)
end

function _order_seed_plan(lattice, kind::Symbol, Q, seed::Integer)
    kind in (:period9, :period27) || throw(ArgumentError("kind must be :period9 or :period27"))
    N = nsites(lattice)
    N % 9 == 0 && lattice.Ly % 3 == 0 ||
        throw(ArgumentError("order seeds require N divisible by 9 and Ly divisible by 3"))
    kind === :period27 && lattice.Lx % 3 != 0 &&
        throw(ArgumentError("period27 seeds require complete three-column motifs (Lx divisible by 3)"))
    Q isa Integer && !(Q isa Bool) && Q == N ÷ 9 ||
        throw(ArgumentError("order seeds implement only the Q=N/9 target sector"))
    !(seed isa Bool) && 0 <= seed <= typemax(Int) ||
        throw(ArgumentError("seed must be a nonnegative integer fitting in Int"))
    # No randomness: seed changes only the motif's lattice origin and the
    # sign of a specified internal pair phase. The two kinds differ in their
    # actual translation symmetries and bond patterns for every origin.
    origin = (mod(Int(seed), 3), mod(fld(Int(seed), 3), 3))
    phase = isodd(seed) ? 0.17 : -0.17
    canonical = kagome_cylinder(lattice.Lx, lattice.Ly)
    lattice.sites == canonical.sites ||
        throw(ArgumentError("order seeds require canonical x,y,A,B,C site ordering"))
    labels = fill("Dn", N)
    pairs = NTuple{2,Int}[]
    spectators, spectator_charges = Int[], Int[]
    expected_sz = zeros(Float64, N)
    for x in 0:(lattice.Lx - 1), y in 0:(lattice.Ly - 1)
        cell = _order_seed_cell(kind, x, y, origin)
        offset = 3 * (x * lattice.Ly + y)
        pair = (offset + cell.pair[1], offset + cell.pair[2])
        spectator = offset + cell.spectator
        labels[pair[1]], labels[pair[2]] = "Up", "Dn"
        labels[spectator] = cell.charge == 1 ? "Up" : "Dn"
        push!(pairs, pair)
        push!(spectators, spectator)
        push!(spectator_charges, cell.charge)
        expected_sz[spectator] = cell.charge / 2
    end
    sum(spectator_charges) == Q || error("motif charge does not match Q")
    count(==("Up"), labels) == (N + Q) ÷ 2 || error("product seed has incorrect up-spin count")
    pair_set = Set(pairs)
    expected_bond = [minmax(b.i, b.j) in pair_set ? -0.25 - 0.5cos(phase) :
                     expected_sz[b.i] * expected_sz[b.j] for b in lattice.bonds]
    translations = kind === :period9 ? [[1, 1], [-1, 2]] : [[3, 0], [0, 3]]
    metadata = Dict{String,Any}(
        "schema_version"=>1, "kind"=>String(kind), "status"=>"unoptimized_density_and_bond_biased_trial",
        "literature_vbc_reproduction"=>false, "hamiltonian_modified"=>false,
        "N"=>N, "Lx"=>lattice.Lx, "Ly"=>lattice.Ly, "Q"=>Int(Q), "physical_M"=>Q/2,
        "wrap_cell_vector"=>[0, lattice.Ly], "site_order"=>"x,y,A,B,C",
        "seed"=>Int(seed), "origin_cell"=>collect(origin), "primitive_cell_translations"=>translations,
        "motif_site_count"=>kind === :period9 ? 9 : 27,
        "seed_rule"=>"origin=(seed mod 3,floor(seed/3) mod 3); phase=+0.17 for odd seed, -0.17 otherwise",
        "spectator_charge_rule"=>"(x-ox-y+oy) mod 3: 0,1 => +1; 2 => -1",
        "pair_rule"=>kind === :period9 ? "(x-ox-y+oy) mod 3 => AB,BC,AC" :
                                        "(x-ox+y-oy) mod 3 => AB,BC,AC",
        "pair_phase_radians"=>phase,
        "pair_state"=>"(|Up_i Dn_j>-exp(i*phi)|Dn_i Up_j>)/sqrt(2), i<j; total physical Sz=0",
        "pair_sites"=>[collect(pair) for pair in pairs], "spectator_sites"=>spectators,
        "spectator_integer_charges"=>spectator_charges,
        "expected_sz"=>expected_sz,
        "bond_endpoint_pairs"=>[[b.i, b.j] for b in lattice.bonds],
        "expected_bond_unweighted_theta0"=>expected_bond,
        "expected_pair_pm_real"=>-0.5cos(phase), "expected_pair_pm_imag"=>-0.5sin(phase),
        "schmidt_rank_bound"=>2, "construction_cutoff"=>0.0,
        "interpretation"=>"Explicit trial biases only; no converged order, ground state, or phase claim")
    return (; labels, pairs, phase, metadata)
end

function _order_seed_pair_gate(si, sj, phase::Real)
    oi, oj, ii, ij = prime(si), prime(sj), dag(si), dag(sj)
    gate = ITensor(ComplexF64, oi, oj, ii, ij)
    coefficient = inv(sqrt(2.0))
    gate[oi=>1, oj=>1, ii=>1, ij=>1] = 1
    gate[oi=>2, oj=>2, ii=>2, ij=>2] = 1
    gate[oi=>1, oj=>2, ii=>1, ij=>2] = coefficient
    gate[oi=>2, oj=>1, ii=>2, ij=>1] = coefficient
    gate[oi=>2, oj=>1, ii=>1, ij=>2] = -cis(phase) * coefficient
    gate[oi=>1, oj=>2, ii=>2, ij=>1] = cis(-phase) * coefficient
    return gate
end

"""
    make_order_seed(lattice, sites; kind=:period9, Q, seed=11)

Return `(psi, metadata)` as a named tuple, with a complex U(1) MPS at exactly
`Q=N/9` and a TOML-compatible, allowlisted metadata dictionary. These are
explicit density/bond-biased trials, not reproductions of named literature
VBCs. No term is added to the Hamiltonian. Use the returned state as `psi0`
in a separately specified NN, theta=0 optimization.

Each primitive three-site cell has one selected intra-cell pair in
`(|Up,Dn>-cis(phi)|Dn,Up>)/sqrt(2)` and a spectator with `2Sz=+1` or `-1`.
The pair has total Sz zero but is not a singlet for the nonzero `phi=±0.17`.
Spectator signs repeat `(+,+,-)` in `(x-y) mod 3`, giving Q=1 per nine sites.

For `:period9`, pair choice AB/BC/AC also depends on `(x-y) mod 3`. Its
primitive translations are `(1,1),(-1,2)` in the `(a1,a2)` basis, with area
three primitive cells (nine sites). For `:period27`, pair choice depends on
`(x+y) mod 3`; jointly the two rules have primitive translations `(3,0),(0,3)`
and area nine primitive cells (27 sites). The finite system is open in x;
these periodicities describe the underlying infinite motif, not an OBC
translation symmetry. Require Ly divisible by three and, for period27, Lx
divisible by three. N27 (3,3) and N54 (6,3) share the same termination.

The seed selects a translated motif and the pair phase sign, without random
noise. Charge-preserving two-spin unitaries prepare the state from a product
MPS. Pairs never leave a cell; the exact Schmidt-rank bound is two, so the
construction's maxdim=2 does not approximate this ansatz. Input indices and
ordering are retained. DMRG must be allowed to grow all needed charge sectors.
"""
function make_order_seed(lattice, sites; kind::Symbol=:period9, Q, seed::Integer=11)
    plan = _order_seed_plan(lattice, kind, Q, seed)
    N = nsites(lattice)
    length(sites) == N || throw(ArgumentError("site count does not match the lattice"))
    all(s -> hasqns(s) && dim(s) == 2 && hastags(s, "S=1/2"), sites) ||
        throw(ArgumentError("expected U(1)-conserving spin-1/2 site indices"))
    all(s -> qn(s=>1) == QN("Sz", 1) && qn(s=>2) == QN("Sz", -1), sites) ||
        throw(ArgumentError("site indices must use integer Sz=2*physical Sz"))
    psi = MPS(ComplexF64, sites, plan.labels)
    gates = [_order_seed_pair_gate(sites[i], sites[j], plan.phase) for (i, j) in plan.pairs]
    psi = apply(gates, psi; cutoff=0.0, maxdim=2, move_sites_back=true)
    normalize!(psi)
    flux(psi) == QN("Sz", Int(Q)) || error("constructed order seed has incorrect charge")
    all(i -> siteind(psi, i) == sites[i], eachindex(sites)) || error("seed construction changed site ordering")
    maxlinkdim(psi) <= 2 || error("constructed seed exceeds its exact Schmidt-rank bound")
    plan.metadata["measured_norm"] = norm(psi)
    plan.metadata["measured_maxlinkdim"] = maxlinkdim(psi)
    return (; psi, metadata=plan.metadata)
end

"""Bounded construction/readout checks; no DMRG, eigensolve, or dense N-spin tensor."""
function order_seed_calibration()
    phase = 0.17
    local_sites = siteinds("S=1/2", 2; conserve_qns=true)
    local_gate = _order_seed_pair_gate(local_sites[1], local_sites[2], phase)
    dense_gate = reshape(Array(local_gate, prime(local_sites[1]), prime(local_sites[2]),
                               dag(local_sites[1]), dag(local_sites[2])), 4, 4)
    unitary_error = norm(dense_gate' * dense_gate - I)
    target = ComplexF64[0, -cis(phase), 1, 0] / sqrt(2)
    pair_error = norm(dense_gate[:, 3] - target)
    @assert max(unitary_error, pair_error) < 2e-14

    # Enumerate the translation stabilizer modulo three independently of
    # finite OBC geometry. Nine-site motif has 3 stabilizers, 27-site only 1.
    stabilizers = Dict{String,Any}()
    for kind in (:period9, :period27)
        shifts = [[dx, dy] for dx in 0:2 for dy in 0:2 if
            all(_order_seed_cell(kind, x, y, (0, 0)) ==
                _order_seed_cell(kind, x + dx, y + dy, (0, 0)) for x in 0:2 for y in 0:2)]
        @assert shifts == (kind === :period9 ? [[0, 0], [1, 1], [2, 2]] : [[0, 0]])
        stabilizers[String(kind)] = shifts
    end
    records = Dict{String,Any}[]
    difference = Dict{String,Any}()
    for Lx in (3, 6)
        lattice = kagome_cylinder(Lx, 3)
        sites = spin_sites(lattice)
        sz_values, bond_values = Vector{Float64}[], Vector{Float64}[]
        for kind in (:period9, :period27)
            result = make_order_seed(lattice, sites; kind, Q=nsites(lattice)÷9, seed=11)
            psi, metadata = result.psi, result.metadata
            measured_sz = sz_profile(psi)
            correlations = spin_correlations(psi)
            measured_bond = [real(correlations.zz[b.i, b.j] + correlations.pm[b.i, b.j])
                             for b in lattice.bonds]
            sz_error = maximum(abs.(measured_sz - metadata["expected_sz"]))
            bond_error = maximum(abs.(measured_bond - metadata["expected_bond_unweighted_theta0"]))
            pm_expected = complex(metadata["expected_pair_pm_real"], metadata["expected_pair_pm_imag"])
            pm_error = maximum(abs(correlations.pm[pair[1], pair[2]] - pm_expected)
                               for pair in metadata["pair_sites"])
            @assert maximum((sz_error, bond_error, pm_error)) < 2e-12
            @assert abs(2sum(measured_sz) - metadata["Q"]) < 2e-12
            @assert abs(norm(psi) - 1) < 2e-12
            @assert abs(imag(pm_expected)) > 0.05 # Physical complex coherence, not just global phase.
            push!(sz_values, measured_sz)
            push!(bond_values, measured_bond)
            push!(records, Dict{String,Any}("N"=>nsites(lattice), "Q"=>metadata["Q"],
                "kind"=>String(kind), "sz_max_error"=>sz_error, "bond_max_error"=>bond_error,
                "complex_pair_coherence_max_error"=>pm_error,
                "norm"=>norm(psi), "maxlinkdim"=>maxlinkdim(psi)))
        end
        sz_difference = maximum(abs.(sz_values[1] - sz_values[2]))
        bond_difference = maximum(abs.(bond_values[1] - bond_values[2]))
        @assert sz_difference > 0.4 && bond_difference > 0.7
        difference[string(nsites(lattice))] = Dict("sz_max_difference"=>sz_difference,
                                                  "bond_max_difference"=>bond_difference)
    end
    return Dict{String,Any}("passed"=>true, "local_gate_unitarity_error"=>unitary_error,
        "local_pair_vector_error"=>pair_error, "translation_stabilizers_mod3"=>stabilizers,
        "measured_cases"=>records, "motif_profile_differences"=>difference,
        "dmrg_run"=>false, "dense_operator_max_dimension"=>4,
        "interpretation"=>"Construction and readout calibration only; no VBC or ground-state claim")
end

if abspath(PROGRAM_FILE) == @__FILE__
    BLAS.set_num_threads(1)
    elapsed = @elapsed metrics = order_seed_calibration()
    println((; julia_version=string(VERSION), julia_threads=Threads.nthreads(),
              blas_threads=BLAS.get_num_threads(), elapsed_seconds=elapsed, metrics))
end
