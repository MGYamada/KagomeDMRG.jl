using Test
using LinearAlgebra
using Random

@testset "Independent sparse eigensolver and exact multiplicity" begin
    for theta in (0.0, 0.37)
        sparse_ref = reference_eigensystem(1, 3, theta; nev=4, blocksize=4)
        dense_ref = eigen(Hermitian(reference_hamiltonian(1, 3, theta).H))
        @test sparse_ref.H isa SparseMatrixCSC
        @test sparse_ref.converged
        @test sparse_ref.solver_info.num_converged >= 4
        @test sparse_ref.values ≈ dense_ref.values[1:4] atol=2e-11 rtol=0
        @test size(sparse_ref.vectors) == (126, 4)
        @test maximum(sparse_ref.residual_norms) <= sparse_ref.settings.tol
        @test sparse_ref.orthogonality_error < 1e-11
        @test sparse_ref.settings.blocksize == 4
        @test sparse_ref.solver_info.algorithm == "KrylovKit.BlockLanczos"
        ground_indices = findall(v -> v - dense_ref.values[1] < 1e-9, dense_ref.values)
        @test length(ground_indices) == (iszero(theta) ? 2 : 1)
        dense_ground = dense_ref.vectors[:, ground_indices]
        sparse_ground = sparse_ref.vectors[:, 1:length(ground_indices)]
        @test norm(sparse_ground * sparse_ground' - dense_ground * dense_ground') < 2e-9
        for j in eachindex(sparse_ref.values)
            @test sparse_ref.residual_norms[j] ≈ norm(sparse_ref.H * sparse_ref.vectors[:, j] -
                sparse_ref.values[j] * sparse_ref.vectors[:, j]) atol=1e-15 rtol=0
        end
    end

    # A fully known spectrum gives an independent multiplicity check, including
    # a complex unitary change of basis. A starting block of size four resolves
    # the three-fold ground eigenspace without relying on roundoff-generated vectors.
    rng = MersenneTwister(81)
    U = Matrix(qr(randn(rng, ComplexF64, 32, 32)).Q)
    spectrum = [-3.0, -3.0, -3.0, -1.0, collect(0.0:27.0)...]
    H = Hermitian(U * Diagonal(spectrum) * U')
    known = _reference_low_eigensystem(H; nev=4, blocksize=4, krylovdim=28)
    @test known.converged
    @test known.values ≈ spectrum[1:4] atol=2e-11 rtol=0
    @test norm(known.vectors[:, 1:3] * known.vectors[:, 1:3]' - U[:, 1:3] * U[:, 1:3]') < 2e-10
    @test maximum(known.residual_norms) <= known.settings.tol

    # Explicit seed changes preserve the complete two-dimensional ground space.
    first_seed = reference_eigensystem(1, 3; seed=45)
    other_seed = reference_eigensystem(1, 3; seed=46, blocksize=4)
    @test first_seed.converged && other_seed.converged
    @test first_seed.values ≈ other_seed.values atol=2e-11 rtol=0
    @test norm(first_seed.vectors[:, 1:2] * first_seed.vectors[:, 1:2]' -
               other_seed.vectors[:, 1:2] * other_seed.vectors[:, 1:2]') < 2e-9

    Random.seed!(923)
    expected_rng_value = rand()
    Random.seed!(923)
    reference_eigensystem(1, 3; nev=1, seed=81)
    @test rand() == expected_rng_value
end

@testset "Sparse eigensolver rejects invalid inputs and reports nonconvergence" begin
    @test_throws ArgumentError reference_eigensystem(0, 3)
    @test_throws ArgumentError reference_eigensystem(1, 2)
    @test_throws ArgumentError reference_eigensystem(3, 3)
    @test_throws ArgumentError reference_eigensystem(1, 4)
    @test_throws ArgumentError reference_eigensystem(typemax(Int), 3)
    @test_throws ArgumentError reference_eigensystem(1, 3, Inf)
    @test_throws ArgumentError reference_eigensystem(1, 3, NaN)
    for kwargs in ((; nev=0), (; nev=126), (; tol=0.0), (; tol=Inf), (; tol=NaN),
                   (; krylovdim=3), (; maxiter=0), (; seed=-1), (; blocksize=0),
                   (; blocksize=64))
        @test_throws ArgumentError reference_eigensystem(1, 3; kwargs...)
    end
    @test_throws ArgumentError _reference_low_eigensystem(zeros(3, 4))
    @test_throws ArgumentError _reference_low_eigensystem(fill(NaN, 8, 8))
    nonhermitian = Matrix{ComplexF64}(I, 8, 8)
    nonhermitian[1, 2] = 1im
    @test_throws ArgumentError _reference_low_eigensystem(nonhermitian)

    failed = reference_eigensystem(1, 3, 0.37; krylovdim=6, maxiter=1, tol=1e-12)
    @test !failed.converged
    @test failed.solver_info.num_converged < failed.settings.nev
    @test failed.solver_info.numiter == 1
    @test length(failed.values) == 3
    @test any(>(failed.settings.tol), failed.residual_norms)
    @test all(isfinite, failed.values)
    @test all(isfinite, failed.residual_norms)
end
