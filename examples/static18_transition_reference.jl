# Independent research-only spin-basis transition reference. No KagomeDMRG,
# production observables, MPO, or Hamiltonian implementation is imported.
# Bounded calibration (largest dense operator: 8 x 8):
# julia --project=. --startup-file=no --threads=1 examples/static18_transition_reference.jl
using LinearAlgebra

function _static18_transition_basis(basis::AbstractVector{<:Integer}, N::Integer)
    !(N isa Bool) && 1 <= N <= 18 ||
        throw(ArgumentError("transition reference requires 1 <= N <= 18"))
    isempty(basis) && throw(ArgumentError("the fixed-charge basis must be nonempty"))
    all(s -> !(s isa Bool) && 0 <= s < (UInt64(1) << N), basis) ||
        throw(ArgumentError("basis states must be nonnegative integers on N sites"))
    states = UInt64.(basis)
    nup = count_ones(first(states))
    all(s -> count_ones(s) == nup, states) ||
        throw(ArgumentError("basis states must belong to one fixed-charge sector"))
    lookup = Dict(state => k for (k, state) in enumerate(states))
    length(lookup) == length(states) || throw(ArgumentError("basis states must be unique"))
    length(states) == binomial(Int(N), nup) ||
        throw(ArgumentError("the fixed-charge basis must be complete"))
    return (; states, lookup, nup)
end

function _static18_transition_vectors(vectors::AbstractMatrix{<:Number}, dimension::Integer)
    size(vectors, 1) == dimension || throw(DimensionMismatch("state and basis lengths differ"))
    size(vectors, 2) > 0 || throw(ArgumentError("at least one state column is required"))
    result = ComplexF64.(vectors)
    all(isfinite, result) || throw(ArgumentError("state amplitudes must be finite"))
    for column in eachcol(result)
        magnitude = norm(column)
        isfinite(magnitude) && magnitude > 0 ||
            throw(ArgumentError("each state must have finite positive norm"))
    end
    return result
end

"""
    transition_amplitudes(ground, excited_matrix, basis, N, bonds)

Return `(sz, bond)` as ComplexF64 matrices of size `N x n_excited` and
`length(bonds) x n_excited`. Entry `[a, alpha]` is `<e_alpha|O_a|g>`.
`excited_matrix` stores one ket per column. The input vectors are preserved;
no normalization or orthogonalization is applied, so these amplitudes remain
sesquilinear in the supplied states. Interpreting sums of squared amplitudes
as spectral weights requires a normalized ground state and orthonormal
excited states, which the caller must independently check.

`O_i=Sz_i`, and each bond uses the *unweighted, zero-flux, isotropic* operator
`Si⋅Sj=Sz_i Sz_j + (S+_i S-_j + S-_i S+_j)/2`. Only a bond's `i,j` fields
are used: stored couplings, family, orientation, and winding are irrelevant
to this definition. In particular this helper is not a finite-flux bond
energy measurement. Basis bit `i-1` is one for Up, as in reference_ed.jl.

Require a complete duplicate-free fixed-Q basis in any order, `N <= 18`,
and finite nonzero states. The implementation applies each operator to a
spin-basis vector and uses Hermitian inner products; no full operator matrix
is constructed. The phase of each transition column depends on the state
phases. For `E_new=E*U`, the amplitude matrix obeys `A_new=A*conj(U)`.
"""
function transition_amplitudes(ground::AbstractVector{<:Number},
        excited_matrix::AbstractMatrix{<:Number}, basis::AbstractVector{<:Integer},
        N::Integer, bonds)
    sector = _static18_transition_basis(basis, N)
    states, lookup = sector.states, sector.lookup
    g = vec(_static18_transition_vectors(reshape(ground, :, 1), length(states)))
    excited = _static18_transition_vectors(excited_matrix, length(states))
    edges = collect(bonds)
    for edge in edges
        hasproperty(edge, :i) && hasproperty(edge, :j) ||
            throw(ArgumentError("bonds must provide i and j fields"))
        all(i -> i isa Integer && !(i isa Bool) && 1 <= i <= N, (edge.i, edge.j)) &&
            edge.i != edge.j || throw(ArgumentError("bond endpoints must be distinct valid integers"))
    end
    count_excited = size(excited, 2)
    sz = zeros(ComplexF64, Int(N), count_excited)
    bond = zeros(ComplexF64, length(edges), count_excited)
    applied = similar(g)
    overlaps = zeros(ComplexF64, count_excited)
    for site in 1:N
        mask = UInt64(1) << (site - 1)
        for (column, state) in enumerate(states)
            applied[column] = (iszero(state & mask) ? -0.5 : 0.5) * g[column]
        end
        mul!(overlaps, excited', applied)
        sz[site, :] = overlaps
    end
    for (edge_index, edge) in enumerate(edges)
        mask_i = UInt64(1) << (edge.i - 1)
        mask_j = UInt64(1) << (edge.j - 1)
        for (column, state) in enumerate(states)
            unlike = iszero(state & mask_i) != iszero(state & mask_j)
            applied[column] = (unlike ? -0.25 : 0.25) * g[column]
            if unlike
                flipped = state ⊻ mask_i ⊻ mask_j
                # The flipped configuration lies in the same fixed-Q sector.
                # Spin 1/2 gives a unit ladder matrix element and coefficient 1/2.
                applied[column] += 0.5 * g[lookup[flipped]]
            end
        end
        mul!(overlaps, excited', applied)
        bond[edge_index, :] = overlaps
    end
    return (; sz, bond)
end

"""
    total_spin_squared(vectors, basis, N) -> Vector{Float64}

Normalized `<S_total^2>` for each supplied state column; a single vector is
also accepted and returns a one-entry vector. In a fixed physical `M=Q/2`
sector, `S² = S⁻S⁺ + Sz(Sz+1)`, hence
`<S²> = M(M+1) + ||S⁺v||²/||v||²`.

Apply `S⁺=sum_i S⁺_i` directly to the adjacent fixed-charge sector, using
unit local raising amplitudes and bit=1 for Up. Work scales with the number
of down spins per configuration, avoiding all-pair exchange sums. The
largest stored object is a state matrix in one fixed-charge sector; no
dense spin-space operator is formed. Normalization uses private copies.
An expectation alone does not establish a definite total-spin eigenstate.
"""
function total_spin_squared(vectors::AbstractMatrix{<:Number},
        basis::AbstractVector{<:Integer}, N::Integer)
    return total_spin_diagnostics(vectors, basis, N).spin_squared
end

function total_spin_squared(vector::AbstractVector{<:Number},
        basis::AbstractVector{<:Integer}, N::Integer)
    return total_spin_squared(reshape(vector, :, 1), basis, N)
end

"""
    total_spin_diagnostics(vectors, basis, N)

Return `(spin_squared, raising_norm_squared)` as two Float64 vectors from
one raising operation. `raising_norm_squared=||S⁺v||²/||v||²` is retained
before adding `M(M+1)`, preserving very small nonzero raising norms that may
be lost to rounding in `<S²>`. The input may be one vector or a column matrix.
Use this combined helper when both diagnostics are needed.
"""
function total_spin_diagnostics(vectors::AbstractMatrix{<:Number},
        basis::AbstractVector{<:Integer}, N::Integer)
    sector = _static18_transition_basis(basis, N)
    states, nup = sector.states, sector.nup
    inputs = _static18_transition_vectors(vectors, length(states))
    for column in eachcol(inputs)
        column ./= norm(column)
    end
    magnetization = nup - Int(N) / 2
    base = magnetization * (magnetization + 1)
    nup == N && return (spin_squared=fill(base, size(inputs, 2)),
                       raising_norm_squared=zeros(Float64, size(inputs, 2)))
    next_states = [state for state in UInt64(0):((UInt64(1) << N) - 1)
                   if count_ones(state) == nup + 1]
    next_lookup = Dict(state => row for (row, state) in enumerate(next_states))
    raised = zeros(ComplexF64, length(next_states), size(inputs, 2))
    for (column, state) in enumerate(states), site in 1:N
        mask = UInt64(1) << (site - 1)
        iszero(state & mask) || continue
        row = next_lookup[state | mask]
        for vector_index in axes(inputs, 2)
            raised[row, vector_index] += inputs[column, vector_index]
        end
    end
    raising_norm_squared = [sum(abs2, view(raised, :, column)) / sum(abs2, view(inputs, :, column))
                            for column in axes(inputs, 2)]
    return (spin_squared=base .+ raising_norm_squared, raising_norm_squared=raising_norm_squared)
end

function total_spin_diagnostics(vector::AbstractVector{<:Number},
        basis::AbstractVector{<:Integer}, N::Integer)
    return total_spin_diagnostics(reshape(vector, :, 1), basis, N)
end

"""Return `||S⁺v||²/||v||²` per state; use total_spin_diagnostics to also obtain `<S²>` in one pass."""
total_spin_raising_norm_squared(vectors, basis, N) =
    total_spin_diagnostics(vectors, basis, N).raising_norm_squared

function _static18_transition_dense_spin(N::Integer, site::Integer, component::Integer)
    1 <= N <= 3 || throw(ArgumentError("dense calibration is limited to three spins"))
    spin = (ComplexF64[0 1; 1 0] / 2,
            ComplexF64[0 -im; im 0] / 2,
            ComplexF64[1 0; 0 -1] / 2)
    identity = Matrix{ComplexF64}(I, 2, 2)
    # Cartesian Pauli basis is (Up, Down), site 1 fastest. Global reference
    # integers instead encode Up with one: calibration embeds by complement.
    return foldl(kron, (i == site ? spin[component] : identity for i in N:-1:1))
end

"""
    transition_reference_calibration() -> Dict

Independent three-spin Cartesian Pauli/Kronecker calibration, including
complex matrix elements, unitary mixing of a two-state target subspace,
projected vectors, `C=conj(A)*transpose(A)`, and `<S²>`. No eigensolve is
performed and dense operators never exceed 8 x 8. Return only TOML-compatible
scalars/arrays, including `passed::Bool`; failed checks throw an assertion.
"""
function transition_reference_calibration()
    N = 3
    basis = UInt64[3, 5, 6]
    omega = cis(2pi / 3)
    ground = ComplexF64[1, omega, omega^2] / sqrt(3)
    excited = hcat(ones(ComplexF64, 3), ComplexF64[1, omega^2, omega]) / sqrt(3)
    saved_ground, saved_excited = copy(ground), copy(excited)
    bonds = [(i=1, j=2, Jxy=2.7, Jz=-1.2, wy=1),
             (i=2, j=3, Jxy=-0.3, Jz=4.1, wy=-1),
             (i=1, j=3, Jxy=0.0, Jz=0.0, wy=0)]
    amplitudes = transition_amplitudes(ground, excited, basis, N, bonds)
    dense_ground = zeros(ComplexF64, 8)
    dense_excited = zeros(ComplexF64, 8, 2)
    for (column, state) in enumerate(basis)
        row = 1 + Int(UInt64(7) ⊻ state)
        dense_ground[row] = ground[column]
        dense_excited[row, :] = excited[column, :]
    end
    spin = [_static18_transition_dense_spin(N, site, component) for site in 1:N, component in 1:3]
    dense_sz = [dot(view(dense_excited, :, alpha), spin[site, 3] * dense_ground)
                for site in 1:N, alpha in 1:2]
    dense_bond = zeros(ComplexF64, length(bonds), 2)
    for (index, edge) in enumerate(bonds)
        operator = sum(spin[edge.i, component] * spin[edge.j, component] for component in 1:3)
        for alpha in 1:2
            dense_bond[index, alpha] = dot(view(dense_excited, :, alpha), operator * dense_ground)
        end
    end
    sz_error = maximum(abs.(amplitudes.sz - dense_sz))
    bond_error = maximum(abs.(amplitudes.bond - dense_bond))
    @assert sz_error < 2e-14 && bond_error < 2e-14
    @assert maximum(abs.(imag.(amplitudes.sz))) > 0.1
    @assert maximum(abs.(imag.(amplitudes.bond))) > 0.1

    cosine, sine, phase = cos(0.43), sin(0.43), cis(0.29)
    mixing = cis(0.17) .* ComplexF64[cosine sine * phase; -sine * conj(phase) cosine]
    @assert norm(mixing' * mixing - I) < 2e-14
    mixed = transition_amplitudes(ground, excited * mixing, basis, N, bonds)
    covariance_error = maximum((maximum(abs.(mixed.sz - amplitudes.sz * conj(mixing))),
                                maximum(abs.(mixed.bond - amplitudes.bond * conj(mixing)))))
    projector_error = 0.0
    weight_error = 0.0
    gram_error = 0.0
    for (original, changed) in ((amplitudes.sz, mixed.sz), (amplitudes.bond, mixed.bond))
        projector_error = max(projector_error,
            maximum(abs.(excited * transpose(original) - (excited * mixing) * transpose(changed))))
        weight_error = max(weight_error, maximum(abs.(sum(abs2, original; dims=2) - sum(abs2, changed; dims=2))))
        gram_error = max(gram_error,
            maximum(abs.(conj(original) * transpose(original) - conj(changed) * transpose(changed))))
    end
    @assert maximum((covariance_error, projector_error, weight_error, gram_error)) < 2e-14

    permutation = [3, 1, 2]
    permuted = transition_amplitudes(ground[permutation], excited[permutation, :],
                                    basis[permutation], N, reverse(bonds))
    @assert maximum(abs.(permuted.sz - amplitudes.sz)) < 2e-14
    @assert maximum(abs.(permuted.bond - reverse(amplitudes.bond; dims=1))) < 2e-14
    scale = 1.7cis(0.38)
    scaled = transition_amplitudes(scale .* ground, excited, basis, N, bonds)
    @assert maximum(abs.(scaled.sz - scale .* amplitudes.sz)) < 2e-14
    @assert maximum(abs.(scaled.bond - scale .* amplitudes.bond)) < 2e-14
    @assert ground == saved_ground && excited == saved_excited

    total_components = [sum(spin[site, component] for site in 1:N) for component in 1:3]
    dense_s2 = sum(component * component for component in total_components)
    # Include a superposition of different S sectors and a nonunit norm state:
    # the helper reports an expectation, rather than inferring S from it.
    vectors = hcat(ground, excited, 1.3 .* ground + 0.4im .* excited[:, 1])
    dense_vectors = hcat(dense_ground, dense_excited,
                         1.3 .* dense_ground + 0.4im .* dense_excited[:, 1])
    expected_s2 = [real(dot(column, dense_s2 * column) / dot(column, column))
                   for column in eachcol(dense_vectors)]
    diagnostics = total_spin_diagnostics(vectors, basis, N)
    measured_s2 = diagnostics.spin_squared
    s2_error = maximum(abs.(measured_s2 - expected_s2))
    dense_raising = total_components[1] + im * total_components[2]
    expected_raising = [sum(abs2, dense_raising * column) / sum(abs2, column)
                        for column in eachcol(dense_vectors)]
    raising_error = maximum(abs.(diagnostics.raising_norm_squared - expected_raising))
    @assert s2_error < 2e-14
    @assert raising_error < 2e-14
    @assert total_spin_raising_norm_squared(vectors, basis, N) ≈ expected_raising atol=2e-14 rtol=0
    @assert measured_s2[1:3] ≈ [0.75, 3.75, 0.75] atol=2e-14 rtol=0
    @assert total_spin_squared(ground, basis, N) ≈ [0.75] atol=2e-14 rtol=0
    @assert total_spin_squared(ComplexF64[2im], UInt64[7], N) ≈ [3.75] atol=2e-14 rtol=0
    @assert total_spin_squared(ComplexF64[2im], UInt64[0], N) ≈ [3.75] atol=2e-14 rtol=0
    @assert ground == saved_ground && excited == saved_excited
    return Dict{String,Any}("passed"=>true, "calibration_N"=>N, "calibration_Q"=>1,
        "dense_operator_dimension"=>8, "reference_basis_dimension"=>length(basis),
        "sz_amplitude_max_error"=>sz_error, "bond_amplitude_max_error"=>bond_error,
        "unitary_covariance_max_error"=>covariance_error,
        "projected_vector_max_error"=>projector_error,
        "projector_gram_max_error"=>gram_error, "subspace_weight_max_error"=>weight_error,
        "total_spin_squared_max_error"=>s2_error,
        "total_spin_raising_norm_squared_max_error"=>raising_error,
        "total_spin_raising_norm_squared_reference"=>expected_raising,
        "total_spin_raising_norm_squared_measured"=>diagnostics.raising_norm_squared,
        "total_spin_squared_reference"=>expected_s2,
        "total_spin_squared_measured"=>measured_s2,
        "basis_permutation_checked"=>true, "raw_amplitude_scaling_checked"=>true,
        "input_preservation_checked"=>true,
        "amplitude_convention"=>"A[a,alpha]=<e_alpha|O_a|g>; no implicit normalization",
        "mixing_convention"=>"E_new=E*U => A_new=A*conj(U)",
        "gram_convention"=>"C=conj(A)*transpose(A)",
        "spin_squared_method"=>"M(M+1)+norm(Splus*v)^2/norm(v)^2")
end

if abspath(PROGRAM_FILE) == @__FILE__
    BLAS.set_num_threads(1)
    elapsed = @elapsed metrics = transition_reference_calibration()
    println((; julia_version=string(VERSION), julia_threads=Threads.nthreads(),
              blas_threads=BLAS.get_num_threads(), elapsed_seconds=elapsed, metrics))
end
