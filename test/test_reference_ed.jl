using Test
using LinearAlgebra
using SparseArrays

_ed_canonical_bond(b) = b.i < b.j ? (b.i, b.j, b.wy, b.Jxy, b.Jz) :
                                                (b.j, b.i, -b.wy, b.Jxy, b.Jz)

@testset "Independent Cartesian-distance bonds" begin
    for (Lx, Ly) in ((1, 3), (2, 3), (3, 4))
        Jxy, Jz = 0.8, 1.2
        expected = reference_bonds(Lx, Ly; Jxy, Jz)
        actual = kagome_cylinder(Lx, Ly; Jxy, Jz).bonds
        @test length(expected) == (6Lx - 2) * Ly
        @test length(unique(_ed_canonical_bond.(expected))) == length(expected)
        @test sort(_ed_canonical_bond.(actual)) == sort(_ed_canonical_bond.(expected))
    end
    @test_throws ArgumentError reference_bonds(1, 2)
end

@testset "Independent spin-basis Hamiltonian" begin
    theta = 0.371
    reference = reference_hamiltonian(1, 3, theta)
    H, basis = reference.H, reference.basis
    @test size(H) == (126, 126)
    @test all(count_ones(s) == 5 for s in basis)
    @test ishermitian(H)
    @test H ≈ reference_hamiltonian(1, 3, theta + 2pi).H atol=2e-14 rtol=0
    @test conj(H) ≈ reference_hamiltonian(1, 3, -theta).H atol=2e-14 rtol=0
    @test norm(imag(H)) > 0.1
    # A(0,0) -> C(0,2) uses the image below the seam: wy = -1.
    # On |down_1,up_9>, the physical S+_1 S-_9 coefficient is exp(-im*theta)/2.
    up_at_2_3_4_5_9 = UInt64(2 + 4 + 8 + 16 + 256)
    up_at_1_2_3_4_5 = UInt64(1 + 2 + 4 + 8 + 16)
    column = findfirst(==(up_at_2_3_4_5_9), basis)
    row = findfirst(==(up_at_1_2_3_4_5), basis)
    @test H[row, column] ≈ cis(-theta) / 2 atol=1e-15 rtol=0
    @test H[column, row] ≈ cis(theta) / 2 atol=1e-15 rtol=0

    full = reference_hamiltonian(1, 3, theta; nup=:all)
    total_q = [2count_ones(s) - 9 for s in full.basis]
    @test norm(Diagonal(total_q) * full.H - full.H * Diagonal(total_q)) < 1e-14
    sector_indices = Int.(basis) .+ 1
    @test full.H[sector_indices, sector_indices] == H

    uniform = reference_hamiltonian(1, 3, theta; gauge=:uniform)
    U = Diagonal(reference_gauge_diagonal(basis, 1, 3, theta))
    @test uniform.basis == basis
    @test uniform.H ≈ U * H * U' atol=2e-14 rtol=0
    @test eigvals(Hermitian(uniform.H)) ≈ eigvals(Hermitian(H)) atol=2e-13 rtol=0

    lat = kagome_cylinder(1, 3)
    chi = gauge_angles(lat, theta)
    reference_U = reference_gauge_diagonal(basis, 1, 3, theta)
    production_U = [cis(sum(chi[i] * _reference_sz_bit(s, i) for i in 1:9)) for s in basis]
    @test production_U ≈ reference_U atol=1e-14 rtol=0
    @test reference_hamiltonian(1, 3, theta; Jxy=0).H ==
          reference_hamiltonian(1, 3, 0; Jxy=0).H
    @test_throws ArgumentError reference_hamiltonian(1, 3; gauge=:invalid)
    @test_throws ArgumentError reference_hamiltonian(1, 3; nup=10)
    @test_throws ArgumentError reference_hamiltonian(2, 3)
end

@testset "Reference observables and sparse construction" begin
    reference = reference_hamiltonian(1, 3, 0.371)
    basis = reference.basis
    product = zeros(ComplexF64, length(basis))
    product[1] = 2im # Check normalization as well as an explicitly complex state.
    @test reference_sz(product, basis, 9) == [fill(0.5, 5); fill(-0.5, 4)]
    @test reference_correlation(product, basis, 9, 1, 9) == -0.25
    @test reference_correlation(product, basis, 9, 1, 1; operators=("S+", "S-")) == 1
    @test reference_correlation(product, basis, 9, 9, 9; operators=("S+", "S-")) == 0
    @test reference_correlation(product, basis, 9, 1, 1; operators=("S-", "S+")) == 0

    eig = eigen(Hermitian(reference.H))
    psi = eig.vectors[:, 1]
    @test sum(reference_sz(psi, basis, 9)) ≈ 0.5 atol=2e-14
    @test reference_correlation(psi, basis, 9, 1, 3; operators=("S+", "S-")) ≈
          conj(reference_correlation(psi, basis, 9, 3, 1; operators=("S+", "S-"))) atol=2e-14
    measured_energy = 0.0 + 0.0im
    for b in reference_bonds(1, 3)
        zz = reference_correlation(psi, basis, 9, b.i, b.j)
        pm = reference_correlation(psi, basis, 9, b.i, b.j; operators=("S+", "S-"))
        measured_energy += b.Jz * zz + b.Jxy * real(cis(b.wy * 0.371) * pm)
    end
    @test measured_energy ≈ eig.values[1] atol=2e-13 rtol=0

    sparse_ref = reference_hamiltonian(2, 3, 0.371; sparse_matrix=true)
    @test sparse_ref.H isa SparseMatrixCSC{ComplexF64, Int}
    @test size(sparse_ref.H) == (43758, 43758)
    @test ishermitian(sparse_ref.H)
    @test all(count_ones(s) == 10 for s in sparse_ref.basis)
end
