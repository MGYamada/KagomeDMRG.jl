# KagomeDMRG local regression tests for the 0.4.31+1 block-rank patch.
@eval module $(gensym())
using LinearAlgebra: Diagonal, Hermitian, norm, svd, eigen, diag
using NDTensors: NDTensors, BlockSparseTensor, array, blockview
using Test: @test, @testset, @test_throws

@testset "Canonical exact-tie selection" begin
    ranks, weights, error = NDTensors._truncate_block_spectrum(
        [[0.5], [0.5]], [(2, 2), (1, 1)]; maxdim = 1, cutoff = 0.0
    )
    @test ranks == [0, 1]
    @test weights == [0.5]
    @test error == 0.5
end

@testset "Per-block constraints and cutoff conventions" begin
    blocks = [[0.8, 0.1], [0.06, 0.04]]
    keys = [(1, 1), (2, 2)]
    @test_throws ArgumentError NDTensors._truncate_block_spectrum(
        blocks, keys; maxdim = 1, min_blockdim = 1, cutoff = 0.0
    )
    cases = (
        (; name="mandatory states exhaust capacity", maxdim=2, cutoff=0.0,
           absolute=false, relative=true, scale=1.0, ranks=[1, 1], error=0.14),
        (; name="one optional slot", maxdim=3, cutoff=0.0,
           absolute=false, relative=true, scale=1.0, ranks=[2, 1], error=0.04),
        (; name="full rank", maxdim=4, cutoff=0.0,
           absolute=false, relative=true, scale=4.0, ranks=[2, 2], error=0.0),
        (; name="relative cumulative budget", maxdim=4, cutoff=0.2,
           absolute=false, relative=true, scale=4.0, ranks=[1, 1], error=0.14),
        (; name="absolute cumulative budget", maxdim=4, cutoff=0.5,
           absolute=false, relative=false, scale=4.0, ranks=[2, 1], error=0.16),
        (; name="absolute individual threshold", maxdim=4, cutoff=0.5,
           absolute=true, relative=true, scale=4.0, ranks=[1, 1], error=0.56))
    @testset "$(case.name)" for case in cases
        scaled = [case.scale .* b for b in blocks]
        ranks, weights, error = NDTensors._truncate_block_spectrum(
            scaled, keys; maxdim=case.maxdim, min_blockdim=1, cutoff=case.cutoff,
            use_absolute_cutoff=case.absolute, use_relative_cutoff=case.relative
        )
        @test ranks == case.ranks
        expected = sort(vcat((scaled[b][1:case.ranks[b]] for b in eachindex(blocks))...); rev=true)
        @test weights == expected
        @test error ≈ case.error atol=1e-14
    end
end

@testset "SVD and eigen apply the same selected block ranks" begin
    cases = ((; name="exact boundary", probabilities=[0.6, 0.2, 0.2, 0.0], minimum_block=0),
             (; name="near boundary", probabilities=[0.6, 0.2001, 0.1999, 0.0], minimum_block=0),
             (; name="mandatory block ranks", probabilities=[0.8, 0.1, 0.06, 0.04], minimum_block=1))
    @testset "$(case.name)" for case in cases
        (; probabilities, minimum_block) = case
        A = BlockSparseTensor{ComplexF64}([(1, 1), (2, 2)], [2, 2], [2, 2])
        R = BlockSparseTensor{ComplexF64}([(1, 1), (2, 2)], [2, 2], [2, 2])
        for b in 1:2
            p = probabilities[(2b - 1):(2b)]
            array(blockview(A, (b, b))) .= Diagonal(sqrt.(p) .* cis.([0.31, -0.47]))
            array(blockview(R, (b, b))) .= Diagonal(p)
        end
        U, S, V, spectrum = svd(A; maxdim = 2, cutoff = 0.0, min_blockdim = minimum_block)
        retained = array(U) * array(S) * transpose(array(V))
        @test size(array(S)) == (2, 2)
        @test length(spectrum.eigs) == 2
        @test norm(array(A) - retained)^2 / norm(array(A))^2 ≈ spectrum.truncerr atol = 1e-14
        @test sum(spectrum.eigs) ≈ norm(retained)^2 atol = 1e-14

        D, W, eigen_spectrum = eigen(Hermitian(R); maxdim = 2, cutoff = 0.0,
            min_blockdim = minimum_block)
        @test size(array(D)) == (2, 2)
        @test length(eigen_spectrum.eigs) == 2
        projected = array(W)' * array(R) * array(W)
        @test sort(real.(diag(projected)); rev = true) ≈ eigen_spectrum.eigs atol = 1e-14
        @test (sum(probabilities) - sum(eigen_spectrum.eigs)) / sum(probabilities) ≈
            eigen_spectrum.truncerr atol = 1e-14
        if minimum_block == 1
            @test_throws ArgumentError svd(A; maxdim = 1, cutoff = 0.0, min_blockdim = 1)
            @test_throws ArgumentError eigen(Hermitian(R); maxdim = 1, cutoff = 0.0, min_blockdim = 1)
        end
    end
end
end
