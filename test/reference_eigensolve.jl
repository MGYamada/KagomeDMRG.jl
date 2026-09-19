# Sparse low-energy reference. Include reference_ed.jl first. The Hamiltonian is
# constructed only by its independent Cartesian-distance/spin-basis code.
using KrylovKit
using Random

function _reference_low_eigensystem(H::AbstractMatrix;
        nev::Integer=3, tol::Real=1e-11, krylovdim::Integer=80,
        maxiter::Integer=200, seed::Integer=1234,
        blocksize::Integer=max(nev, 2))
    dimension = size(H, 1)
    size(H, 2) == dimension || throw(ArgumentError("H must be square"))
    2 <= dimension <= 43758 || throw(ArgumentError("reference dimension must lie in 2:43758"))
    1 <= nev < dimension || throw(ArgumentError("nev must lie in 1:(dimension-1)"))
    1 <= blocksize <= dimension ÷ 2 ||
        throw(ArgumentError("blocksize must lie in 1:div(dimension,2)"))
    isfinite(tol) && tol > 0 || throw(ArgumentError("tol must be finite and positive"))
    maxiter >= 1 || throw(ArgumentError("maxiter must be positive"))
    seed >= 0 || throw(ArgumentError("seed must be nonnegative"))
    effective_krylovdim = min(krylovdim, dimension)
    effective_krylovdim >= max(2blocksize, nev + blocksize) ||
        throw(ArgumentError("krylovdim must allow at least two blocks and nev+blocksize vectors"))
    all(isfinite, H isa SparseMatrixCSC ? nonzeros(H) : H) ||
        throw(ArgumentError("H entries must be finite"))
    ishermitian(H) || throw(ArgumentError("H must be Hermitian"))
    settings = (; nev=Int(nev), tol=Float64(tol),
        krylovdim=Int(effective_krylovdim), maxiter=Int(maxiter), seed=Int(seed),
        blocksize=Int(blocksize), orthogonality_tolerance=sqrt(eps(Float64)))
    rng = MersenneTwister(seed)
    initial = KrylovKit.Block([randn(rng, ComplexF64, dimension) for _ in 1:blocksize])
    algorithm = KrylovKit.BlockLanczos(; tol=settings.tol,
        krylovdim=settings.krylovdim, maxiter=settings.maxiter,
        qr_tol=eps(Float64), verbosity=0)
    all_values, all_vectors, info = KrylovKit.eigsolve(H, initial, Int(nev), :SR, algorithm)
    selection = sortperm(all_values)[1:min(nev, length(all_values))]
    values = Float64.(all_values[selection])
    vectors = hcat(all_vectors[selection]...)
    # Recompute with the full independent matrix instead of accepting the
    # solver's projected residual estimate as a certificate.
    residual_norms = [norm(H * vectors[:, j] - values[j] * vectors[:, j])
                      for j in eachindex(values)]
    orthogonality_error = norm(vectors' * vectors - I)
    converged = length(values) == nev && info.converged >= nev &&
                all(isfinite, residual_norms) && all(<=(settings.tol), residual_norms) &&
                orthogonality_error <= settings.orthogonality_tolerance
    solver_info = (; algorithm="KrylovKit.BlockLanczos",
        version=string(pkgversion(KrylovKit)), num_converged=info.converged,
        numiter=info.numiter, numops=info.numops,
        returned_count=length(all_values),
        reported_residual_norms=Float64.(info.normres[selection]))
    return (; values, vectors, residual_norms, converged, settings,
        solver_info, orthogonality_error)
end

"""
    reference_eigensystem(Lx, Ly, theta=0.0;
        nev=3, tol=1e-11, krylovdim=80, maxiter=200, seed=1234,
        blocksize=max(nev, 2))

Find the lowest `nev` eigenpairs of the independent, sparse, nearest-neighbor
isotropic Heisenberg Hamiltonian in the `M/Msat=1/9` sector, in seam gauge.
At most 18 sites (dimension 43758) are supported. Include `reference_ed.jl`
before this file; no production MPO or DMRG routine is used.

The returned named tuple contains `H`, `basis`, sorted `values`, a column matrix
`vectors`, independently computed `residual_norms`, `orthogonality_error`,
`converged::Bool`, `settings`, and `solver_info`. A failed eigensolve returns
its Ritz pairs with `converged=false`; callers must check this flag before
using them as a reference. Convergence requires both KrylovKit's convergence
count and each full-matrix residual to meet the requested absolute `tol`.

The complex random starting block uses a private RNG. Unlike single-vector
Lanczos, block Lanczos can resolve exact multiplicity up to `blocksize`.
This is still a partial spectrum: residual convergence does not prove that no
lower eigenvector was missed, or that the returned ground space is complete.
Compare another seed/block size and inspect the last requested level before
claiming multiplicity. A finite-cluster level spacing is not a bulk gap.

The block API and limitation follow KrylovKit 0.10.4's installed source and
https://jutho.github.io/KrylovKit.jl/stable/man/eig/ .
"""
function reference_eigensystem(Lx::Integer, Ly::Integer, theta::Real=0.0;
        nev::Integer=3, tol::Real=1e-11, krylovdim::Integer=80,
        maxiter::Integer=200, seed::Integer=1234,
        blocksize::Integer=max(nev, 2))
    Lx >= 1 && Ly >= 3 || throw(ArgumentError("require Lx >= 1 and Ly >= 3"))
    # Bound before multiplying integers or enumerating spin states.
    Ly <= 6 && Lx <= 6 ÷ Ly || throw(ArgumentError("the reference supports at most 18 sites"))
    (3Lx * Ly) % 9 == 0 || throw(ArgumentError("the 1/9 target requires N divisible by 9"))
    isfinite(theta) || throw(ArgumentError("theta must be finite"))
    reference = reference_hamiltonian(Lx, Ly, theta; sparse_matrix=true)
    eigensystem = _reference_low_eigensystem(reference.H;
        nev, tol, krylovdim, maxiter, seed, blocksize)
    return (; reference..., eigensystem...)
end
