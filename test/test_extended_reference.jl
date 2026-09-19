@testset "Extended exchange versus planar geometry" begin
    @testset "Lx=$Lx Ly=$Ly" for (Lx, Ly) in ((1, 4), (3, 4))
        lattice = kagome_j1j2j3_cylinder(Lx, Ly; J1=0.8, J2=0.3, J3=0.6)
        reference = reference_extended_bonds(Lx, Ly; J1=0.8, J2=0.3, J3=0.6)
        actual = [b.i < b.j ? (b.i, b.j, b.wy, b.Jxy, b.Jz, family) :
                              (b.j, b.i, -b.wy, b.Jxy, b.Jz, family)
                  for (b, family) in zip(lattice.bonds, bond_families(lattice))]
        expected = [(b.i, b.j, b.wy, b.Jxy, b.Jz, b.family) for b in reference]
        @test sort(actual) == sort(expected)
        @test Tuple(count(b -> b.family === family, reference) for family in (:J1, :J2, :J3)) ==
              ((6Lx - 2)Ly, (6Lx - 4)Ly, (3Lx - 2)Ly)
        if Lx == 1
            # No hexagon fits inside a single axial cell. Its valid opposite
            # endpoints nevertheless yield Ly J3 bonds at the open ends.
            plane = _reference_extended_hexagons(Lx, Ly)
            @test !any(all(0 <= plane.vertices[i].x < Lx for i in cycle)
                       for cycle in plane.hexagons)
        else
            # A(0,0)-B(0,0)-A(1,0) is a straight NN chain: its endpoints
            # are distance one, but must not be classified as hexagon J3.
            a, middle, b = 1, 2, 3Ly + 1
            nn = Set((e.i, e.j, e.wy) for e in reference if e.family === :J1)
            j3 = Set((e.i, e.j, e.wy) for e in reference if e.family === :J3)
            @test (a, middle, 0) in nn && (middle, b, 0) in nn
            @test !((a, b, 0) in j3)
        end
    end
end

@testset "Complex extended Hamiltonian reference" begin
    theta = 0.371
    couplings = (; J1=0.8, J2=0.3, J3=0.6)
    reference = reference_extended_hamiltonian(1, 3, theta; Q=1, couplings...)
    @test ishermitian(reference.H)
    @test norm(imag(reference.H)) > 0.1
    periodic = reference_extended_hamiltonian(1, 3, theta + 2pi; Q=1, couplings...)
    @test reference.H ≈ periodic.H atol=2e-14 rtol=0
    uniform = reference_extended_hamiltonian(1, 3, theta; Q=1, gauge=:uniform, couplings...)
    U = Diagonal(reference_gauge_diagonal(reference.basis, 1, 3, theta))
    @test uniform.H ≈ U * reference.H * U' atol=2e-14 rtol=0

    lattice = kagome_j1j2j3_cylinder(1, 3; couplings...)
    sites = spin_sites(lattice)
    full = full_mpo_matrix(twisted_exchange_mpo(sites, lattice, theta), sites)
    indices = 2^9 .- Int.(reference.basis)
    @test full[indices, indices] ≈ reference.H atol=2e-12 rtol=0
    nn_only = reference_extended_hamiltonian(1, 3, theta; J2=0, J3=0)
    @test nn_only.H == reference_hamiltonian(1, 3, theta).H
    @test_throws ArgumentError reference_extended_hamiltonian(2, 3; Q=0)
    @test_throws ArgumentError reference_extended_hamiltonian(typemax(Int), 3; Q=0)
end
