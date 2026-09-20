# Research-only spin-basis scalar-chirality reference. This file imports no
# KagomeDMRG, ITensor, MPO, or production ladder-operator implementation.
# Run the bounded analytic checks with:
# julia --project=. --startup-file=no --threads=1 examples/extended18_chirality_reference.jl
using LinearAlgebra

function _extended18_cartesian_chirality()
    spin = (ComplexF64[0 1; 1 0] / 2,
            ComplexF64[0 -im; im 0] / 2,
            ComplexF64[1 0; 0 -1] / 2)
    operator = zeros(ComplexF64, 8, 8)
    for a in 1:3, b in 1:3, c in 1:3
        epsilon = (a - b) * (b - c) * (c - a) / 2
        # Vertex 1 is the fastest index; each Pauli basis is (Up, Down).
        operator .+= epsilon .* kron(spin[c], spin[b], spin[a])
    end
    return operator
end

"""
    extended18_reference_chiralities(v, basis, N, triangles, theta)

Normalized seam-gauge scalar chirality in the iteration order of `triangles`.
Each triangle supplies ordered `sites` and integer periodic `image_y` fields;
the supplied order determines the orientation. For CCW triangles the operator
is `Si⋅(Sj×Sk)`, locally dressed by `U*C*U†`, where
`U=exp(im*sum(alpha[a]*Sz_a))` and `alpha[a]=-image_y[a]*theta`.

The independent local operator uses Cartesian Pauli matrices and the
Levi-Civita tensor. Its local bit is one for Down, whereas `basis` follows
`test/reference_ed.jl`: bit `i-1` is one for Up at physical site `i`.
Embedding explicitly translates between these conventions and uses a bit
lookup; it never constructs an operator larger than 8×8 or a dense N-spin
tensor. Images with the same tuple share a local matrix within this call.

Require `3 <= N <= 18`, a complete, duplicate-free fixed-charge basis in any
order, a matching finite nonzero vector, and finite real `theta`. Values are
computed in ComplexF64 and returned as Float64. The input is not modified.
This helper does not validate geometric CCW orientation or perform an
eigensolve; its caller must independently validate the supplied triangles
and ground-state accuracy. Finite chirality does not identify a CSL phase.
"""
function extended18_reference_chiralities(v::AbstractVector{<:Number},
        basis::AbstractVector{<:Integer}, N::Integer, triangles, theta::Real)
    !(N isa Bool) && 3 <= N <= 18 ||
        throw(ArgumentError("the chirality reference requires 3 <= N <= 18"))
    length(v) == length(basis) || throw(DimensionMismatch("state and basis lengths differ"))
    isempty(basis) && throw(ArgumentError("the fixed-charge basis must be nonempty"))
    angle = Float64(theta)
    isfinite(angle) || throw(ArgumentError("theta must be finite"))
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
    state_vector = ComplexF64.(v)
    all(isfinite, state_vector) || throw(ArgumentError("state amplitudes must be finite"))
    state_norm = norm(state_vector)
    isfinite(state_norm) && state_norm > 0 ||
        throw(ArgumentError("state must have finite positive norm"))
    state_vector ./= state_norm
    norm_squared = sum(abs2, state_vector)

    bare = _extended18_cartesian_chirality()
    cache = Dict{NTuple{3,Int},Matrix{ComplexF64}}()
    values = Float64[]
    for triangle in triangles
        hasproperty(triangle, :sites) && hasproperty(triangle, :image_y) ||
            throw(ArgumentError("triangles must supply sites and image_y fields"))
        length(triangle.sites) == length(triangle.image_y) == 3 ||
            throw(ArgumentError("each triangle must have three sites and images"))
        all(i -> i isa Integer && !(i isa Bool) && 1 <= i <= N, triangle.sites) &&
            length(unique(triangle.sites)) == 3 ||
            throw(ArgumentError("triangle sites must be distinct valid integers"))
        all(m -> m isa Integer && !(m isa Bool) && typemin(Int) <= m <= typemax(Int),
            triangle.image_y) || throw(ArgumentError("triangle images must fit in Int"))
        vertices = ntuple(a -> Int(triangle.sites[a]), 3)
        images = ntuple(a -> Int(triangle.image_y[a]), 3)
        operator = get!(cache, images) do
            alpha = ntuple(a -> -Float64(images[a]) * angle, 3)
            all(isfinite, alpha) || throw(ArgumentError("triangle angles must be finite"))
            phases = [sum(alpha[a] * (isodd(bits >> (a - 1)) ? -0.5 : 0.5)
                          for a in 1:3) for bits in 0:7]
            all(isfinite, phases) || throw(ArgumentError("local rotation phases must be finite"))
            rotation = Diagonal(cis.(phases))
            rotation * bare * rotation'
        end
        masks = ntuple(a -> UInt64(1) << (vertices[a] - 1), 3)
        triangle_mask = masks[1] | masks[2] | masks[3]
        # Row bits index Down spins locally; the embedded state uses Up bits.
        lifted = [foldl(|, (isodd(bits >> (a - 1)) ? UInt64(0) : masks[a]
                           for a in 1:3); init=UInt64(0)) for bits in 0:7]
        expectation = 0.0 + 0.0im
        for (column, state) in enumerate(states)
            amplitude = state_vector[column]
            iszero(amplitude) && continue
            local_column = 1 + sum(iszero(state & masks[a]) ? (1 << (a - 1)) : 0
                                   for a in 1:3)
            spectators = state & ~triangle_mask
            for local_row in 1:8
                element = operator[local_row, local_column]
                iszero(element) && continue
                # Cartesian construction conserves charge exactly. A missing
                # state here would expose an embedding/operator inconsistency.
                row = lookup[spectators | lifted[local_row]]
                expectation += conj(state_vector[row]) * element * amplitude
            end
        end
        expectation /= norm_squared
        isfinite(expectation) && abs(imag(expectation)) <= 1e-10 * max(1, abs(real(expectation))) ||
            error("chirality expectation is nonfinite or has a significant imaginary part")
        push!(values, real(expectation))
    end
    return values
end

"""Small analytic checks only: no DMRG, Hamiltonian construction, or eigensolve."""
function extended18_reference_chirality_selfcheck()
    bare = _extended18_cartesian_chirality()
    charge = Diagonal([3 - 2count_ones(bits) for bits in 0:7])
    @assert norm(bare - bare') < 2e-14
    @assert norm(charge * bare - bare * charge) < 2e-14

    omega = cis(2pi / 3)
    basis = UInt64[3, 5, 6] # Down at vertex 3, 2, 1, respectively.
    v = ComplexF64[omega^2, omega, 1] / sqrt(3)
    triangle = (sites=(1, 2, 3), image_y=(0, 0, 0))
    expected = sqrt(3) / 4
    positive = only(extended18_reference_chiralities(v, basis, 3, [triangle], 0.0))
    negative = only(extended18_reference_chiralities(conj(v), basis, 3, [triangle], 0.0))
    @assert abs(positive - expected) < 2e-14
    @assert abs(negative + expected) < 2e-14
    cyclic = (sites=(2, 3, 1), image_y=(0, 0, 0))
    reversed = (sites=(1, 3, 2), image_y=(0, 0, 0))
    @assert maximum(abs.(extended18_reference_chiralities(v, basis, 3, [cyclic, reversed], 0.0) -
                         [expected, -expected])) < 2e-14
    @assert only(extended18_reference_chiralities([1.0], UInt64[7], 3, [triangle], 0.0)) == 0

    theta = 0.37
    dressed = (sites=(1, 2, 3), image_y=(1, -2, 3))
    alpha = -theta .* dressed.image_y
    rotation = [cis(sum(alpha[a] * (isodd(state >> (a - 1)) ? 0.5 : -0.5)
                         for a in 1:3)) for state in basis]
    rotated = rotation .* v
    before = copy(rotated)
    generic = only(extended18_reference_chiralities(2cis(0.19) .* rotated, basis, 3,
                                                  [dressed], theta))
    @assert abs(generic - expected) < 2e-14
    flipped = (sites=(1, 3, 2), image_y=(1, 3, -2))
    shifted = (sites=(1, 2, 3), image_y=(5, 2, 7))
    @assert maximum(abs.(extended18_reference_chiralities(rotated, basis, 3,
                         [flipped, shifted], theta) - [-expected, expected])) < 2e-14
    @assert only(extended18_reference_chiralities(reverse(rotated), reverse(basis), 3,
                                                [dressed], theta)) ≈ expected
    @assert rotated == before

    # A three-component chiral state embedded at a nonlocal seam triangle in
    # the complete N18,Q0 sector checks indexing without a large dense tensor.
    basis18 = [s for s in UInt64(0):((UInt64(1) << 18) - 1) if count_ones(s) == 9]
    lookup18 = Dict(s => i for (i, s) in enumerate(basis18))
    seam = (sites=(10, 2, 18), image_y=(0, 0, -1))
    spectator_up = setdiff(1:18, collect(seam.sites))[1:7]
    base = foldl(|, (UInt64(1) << (i - 1) for i in spectator_up); init=UInt64(0))
    amplitudes = ComplexF64[1, omega, omega^2] / sqrt(3)
    embedded_errors = Float64[]
    for flux in (0.0, theta)
        v18 = zeros(ComplexF64, length(basis18))
        for down in 1:3
            state = base | foldl(|, (a == down ? UInt64(0) : UInt64(1) << (seam.sites[a] - 1)
                                   for a in 1:3); init=UInt64(0))
            phase = sum(-seam.image_y[a] * flux * (a == down ? -0.5 : 0.5) for a in 1:3)
            v18[lookup18[state]] = amplitudes[down] * cis(phase)
        end
        measured = only(extended18_reference_chiralities(v18, basis18, 18, [seam], flux))
        push!(embedded_errors, abs(measured - expected))
    end
    @assert maximum(embedded_errors) < 2e-14

    function rejects(f, expected_type)
        try
            f()
        catch error
            return error isa expected_type
        end
        return false
    end
    @assert rejects(() -> extended18_reference_chiralities(zero(v), basis, 3, [triangle], 0.0), ArgumentError)
    @assert rejects(() -> extended18_reference_chiralities([NaN, 0, 0], basis, 3, [triangle], 0.0), ArgumentError)
    @assert rejects(() -> extended18_reference_chiralities(v, UInt64[3, 5, 5], 3, [triangle], 0.0), ArgumentError)
    @assert rejects(() -> extended18_reference_chiralities(v[1:2], basis[1:2], 3, [triangle], 0.0), ArgumentError)
    @assert rejects(() -> extended18_reference_chiralities(v, UInt64[1, 5, 6], 3, [triangle], 0.0), ArgumentError)
    @assert rejects(() -> extended18_reference_chiralities(v, basis, 3, [triangle], NaN), ArgumentError)
    @assert rejects(() -> extended18_reference_chiralities(v[1:2], basis, 3, [triangle], 0.0), DimensionMismatch)
    return (; positive, negative, generic_angle_error=abs(generic - expected),
            embedded18_max_error=maximum(embedded_errors), basis18_dimension=length(basis18))
end

if abspath(PROGRAM_FILE) == @__FILE__
    BLAS.set_num_threads(1)
    elapsed = @elapsed metrics = extended18_reference_chirality_selfcheck()
    println((; status="passed", julia_version=string(VERSION), julia_threads=Threads.nthreads(),
              blas_threads=BLAS.get_num_threads(), elapsed_seconds=elapsed, metrics))
end
