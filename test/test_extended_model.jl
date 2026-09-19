@testset "Hexagon J1-J2-J3 cylinder geometry" begin
    for (lx, ly) in ((1, 4), (2, 3))
        lattice = kagome_j1j2j3_cylinder(lx, ly; J1=1.3, J2=-0.2, J3=0.7)
        families = bond_families(lattice)
        counts = ((6lx - 2) * ly, (6lx - 4) * ly, (3lx - 2) * ly)
        for (family, number, coupling, winding) in zip(
                (:J1, :J2, :J3), counts, (1.3, -0.2, 0.7),
                (2lx - 1, 4lx - 2, 2lx - 1))
            bonds = lattice.bonds[families .== family]
            @test length(bonds) == number
            @test all(b -> b.Jxy == coupling && b.Jz == coupling, bonds)
            @test count(b -> b.wy == -1, bonds) == winding
        end
        @test all(b -> b.wy in (-1, 0), lattice.bonds)
        images = Set(b.i < b.j ? (b.i, b.j, b.wy) : (b.j, b.i, -b.wy)
                     for b in lattice.bonds)
        @test length(images) == sum(counts)
        reversed = KagomeCylinder(lx, ly, lattice.sites, reverse_bond.(lattice.bonds))
        @test bond_families(reversed) == families

        reduced = kagome_j1j2j3_cylinder(lx, ly; J1=0.7, J2=0, J3=0)
        nn = kagome_cylinder(lx, ly; Jxy=0.7, Jz=0.7)
        @test filter(b -> !iszero(b.Jxy), reduced.bonds) == nn.bonds
        @test length(reduced.bonds) == sum(counts)
        @test bond_families(reduced) == families
    end

    lattice = kagome_j1j2j3_cylinder(2, 3)
    b = first(lattice.bonds)
    straight_chain = Bond(site_index(lattice, 0, 0, :A),
                          site_index(lattice, 1, 0, :A), 0.5, 0.5, 0)
    probes = [Bond(b.i, b.j, 7.0, -2.0, b.wy), straight_chain,
              Bond(b.i, b.j, b.Jxy, b.Jz, b.wy + 1)]
    edited = KagomeCylinder(lattice.Lx, lattice.Ly, lattice.sites, probes)
    @test bond_families(edited) == [:J1, :other, :other]

    for coupling in ((; J1=NaN), (; J2=Inf), (; J3=-Inf))
        @test_throws ArgumentError kagome_j1j2j3_cylinder(1, 4; coupling...)
    end
    @test_throws ArgumentError kagome_j1j2j3_cylinder(0, 4)
    @test_throws ArgumentError kagome_j1j2j3_cylinder(1, 2)
end
