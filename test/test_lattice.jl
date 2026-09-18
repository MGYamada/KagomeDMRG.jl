function lattice_oriented_edge(lattice, i, j)
    for b in lattice.bonds
        b.i == i && b.j == j && return b
        b.i == j && b.j == i && return reverse_bond(b)
    end
    error("test loop traverses a missing edge: $i -> $j")
end

function lattice_loop_edges(lattice, vertices)
    return [lattice_oriented_edge(lattice, vertices[k], vertices[mod1(k + 1, length(vertices))])
            for k in eachindex(vertices)]
end

@testset "Kagome geometry and site ordering" begin
    for lx in (1, 2, 4), ly in (3, 4, 6)
        lattice = kagome_cylinder(lx, ly)
        @test nsites(lattice) == 3 * lx * ly
        @test length(lattice.bonds) == (6 * lx - 2) * ly
        @test length(Set((min(b.i, b.j), max(b.i, b.j)) for b in lattice.bonds)) ==
              length(lattice.bonds)
        @test all(b -> b.Jxy == 1.0 && b.Jz == 1.0, lattice.bonds)
        @test count(b -> b.wy == -1, lattice.bonds) == 2 * lx - 1
        @test all(b -> b.wy in (-1, 0), lattice.bonds)

        degrees = zeros(Int, nsites(lattice))
        wrap = (ly / 2, ly * sqrt(3.0) / 2)
        for b in lattice.bonds
            degrees[b.i] += 1
            degrees[b.j] += 1
            ri, rj = lattice.sites[b.i].position, lattice.sites[b.j].position
            displacement = ntuple(d -> rj[d] + b.wy * wrap[d] - ri[d], 2)
            @test sum(abs2, displacement) ≈ 1 / 4 atol=2e-14
        end
        # This edge termination removes one neighbor from A,C on the left
        # and two neighbors from B on the right (including the Lx=1 case).
        for s in lattice.sites
            expected_degree = s.sublattice === :B ? 2 + 2 * (s.x < lx - 1) :
                              3 + (s.x > 0)
            @test degrees[s.index] == expected_degree
            @test s.index == 3 * (s.x * ly + s.y) +
                  (s.sublattice === :A ? 1 : s.sublattice === :B ? 2 : 3)
            @test site_index(lattice, s.x, s.y, s.sublattice) == s.index
        end

        # Independent geometric oracle: every pair at nearest-neighbor
        # distance in any adjacent periodic image appears exactly once,
        # with a winding consistent with that image.
        expected_images = Set{Tuple{Int,Int,Int}}()
        for i in 1:nsites(lattice), j in (i + 1):nsites(lattice), wy in -1:1
            ri, rj = lattice.sites[i].position, lattice.sites[j].position
            distance2 = sum((rj[d] + wy * wrap[d] - ri[d])^2 for d in 1:2)
            isapprox(distance2, 1 / 4; atol=2e-14) && push!(expected_images, (i, j, wy))
        end
        actual_images = Set(b.i < b.j ? (b.i, b.j, b.wy) : (b.j, b.i, -b.wy)
                            for b in lattice.bonds)
        @test actual_images == expected_images
    end

    anisotropic = kagome_cylinder(2, 3; Jxy=0.7, Jz=1.2)
    @test all(b -> b.Jxy == 0.7 && b.Jz == 1.2, anisotropic.bonds)
    @test site_index(anisotropic, 0, 0, :A) == 1
    @test site_index(anisotropic, 1, 2, :C) == 18
    @test anisotropic.sites[1].position == (0.0, 0.0)
    @test anisotropic.sites[2].position == (0.5, 0.0)
    @test anisotropic.sites[3].position == (0.25, sqrt(3.0) / 4)
    @test_throws ArgumentError kagome_cylinder(0, 3)
    @test_throws ArgumentError kagome_cylinder(-1, 3)
    @test_throws ArgumentError kagome_cylinder(1, 2)
    @test_throws ArgumentError kagome_cylinder(1, 1)
    @test_throws ArgumentError kagome_cylinder(1, 0)
    @test_throws ArgumentError kagome_cylinder(1, 3; Jxy=Inf)
    @test_throws ArgumentError kagome_cylinder(1, 3; Jz=NaN)
    @test_throws ArgumentError site_index(anisotropic, -1, 0, :A)
    @test_throws ArgumentError site_index(anisotropic, 2, 0, :A)
    @test_throws ArgumentError site_index(anisotropic, 0, -1, :A)
    @test_throws ArgumentError site_index(anisotropic, 0, 3, :A)
    @test_throws ArgumentError site_index(anisotropic, 0, 0, :D)
end

@testset "Target charge and geometric cuts" begin
    for n in (9, 18, 27, 108, 324)
        sector = target_sector(n)
        @test sector.N == n
        @test sector.Q == n ÷ 9
        @test sector.M == n / 18
        @test sector.Nup + sector.Ndown == n
        @test sector.Nup - sector.Ndown == sector.Q == 2 * sector.M
        @test sector.M / (n / 2) ≈ 1 / 9
    end
    @test target_sector(9).M == 0.5
    @test target_sector(18).M == 1.0
    for n in (-9, 0, 1, 8, 10, 17)
        @test_throws ArgumentError target_sector(n)
    end
    lattice = kagome_cylinder(4, 3)
    for cut in 1:3
        indices = right_region(lattice, cut)
        @test indices == collect((9 * cut + 1):36)
        @test length(indices) == 9 * (4 - cut)
        @test all(i -> lattice.sites[i].x >= cut, indices)
        @test all(i -> lattice.sites[i].x < cut, setdiff(1:36, indices))
    end
    @test_throws ArgumentError right_region(lattice, 0)
    @test_throws ArgumentError right_region(lattice, 4)
    @test_throws ArgumentError right_region(kagome_cylinder(1, 3), 1)
end

@testset "Signed winding and contractible loops" begin
    lattice = kagome_cylinder(3, 4)
    index(x, y, s) = site_index(lattice, x, mod(y, lattice.Ly), s)
    loops = Vector{Int}[]
    for x in 0:2, y in 0:3
        push!(loops, [index(x, y, :A), index(x, y, :B), index(x, y, :C)])
        if x > 0
            push!(loops, [index(x, y, :A), index(x - 1, y, :B), index(x, y - 1, :C)])
        end
        if x < 2
            push!(loops, [index(x, y, :A), index(x, y, :B),
                          index(x + 1, y - 1, :C), index(x + 1, y - 1, :A),
                          index(x, y - 1, :B), index(x, y - 1, :C)])
        end
    end
    @test count(v -> length(v) == 3, loops) == 20
    @test count(v -> length(v) == 6, loops) == 8
    theta = 0.731
    for vertices in loops
        edges = lattice_loop_edges(lattice, vertices)
        @test sum(b.wy for b in edges) == 0
        for gauge in (:seam, :uniform)
            @test sum(bond_phase(lattice, b, theta; gauge) for b in edges) ≈ 0 atol=1e-14
        end
    end
    for x in 0:2
        circumference = [index(x, y, s) for y in 0:3 for s in (:A, :C)]
        edges = lattice_loop_edges(lattice, circumference)
        backwards = lattice_loop_edges(lattice, reverse(circumference))
        @test sum(b.wy for b in edges) == 1
        @test sum(b.wy for b in backwards) == -1
        for gauge in (:seam, :uniform)
            @test sum(bond_phase(lattice, b, theta; gauge) for b in edges) ≈ theta
            @test sum(bond_phase(lattice, b, theta; gauge) for b in backwards) ≈ -theta
        end
    end

    for b in lattice.bonds
        reversed = reverse_bond(b)
        @test (reversed.i, reversed.j, reversed.wy) == (b.j, b.i, -b.wy)
        @test (reversed.Jxy, reversed.Jz) == (b.Jxy, b.Jz)
        @test reverse_bond(reversed) == b
        for gauge in (:seam, :uniform)
            angle = bond_phase(lattice, b, theta; gauge)
            @test bond_phase(lattice, reversed, theta; gauge) ≈ -angle
        end
        @test cis(bond_phase(lattice, b, theta + 2pi)) ≈ cis(bond_phase(lattice, b, theta))
    end

    chi = gauge_angles(lattice, theta)
    @test length(chi) == nsites(lattice)
    @test chi[index(0, 0, :A)] == 0.0
    @test chi[index(0, 0, :C)] ≈ -theta / (2 * lattice.Ly)
    for b in lattice.bonds
        @test bond_phase(lattice, b, theta; gauge=:uniform) ≈
              bond_phase(lattice, b, theta) + chi[b.i] - chi[b.j] atol=1e-15
        # A uniform gauge depends only on the local unwrapped a2 displacement.
        si, sj = lattice.sites[b.i], lattice.sites[b.j]
        dy = sj.y - si.y + b.wy * lattice.Ly +
             (sj.sublattice === :C ? 0.5 : 0.0) - (si.sublattice === :C ? 0.5 : 0.0)
        @test bond_phase(lattice, b, theta; gauge=:uniform) ≈ theta * dy / lattice.Ly atol=1e-15
    end
    @test_throws ArgumentError Bond(0, 2, 1.0, 1.0, 0)
    @test_throws ArgumentError Bond(1, -2, 1.0, 1.0, 0)
    @test_throws ArgumentError Bond(1, 1, 1.0, 1.0, 1)
    @test_throws ArgumentError Bond(1, 2, NaN, 1.0, 0)
    @test_throws ArgumentError Bond(1, 2, 1.0, Inf, 0)
    @test_throws ArgumentError bond_phase(lattice, first(lattice.bonds), theta; gauge=:unknown)
    @test_throws ArgumentError bond_phase(lattice, Bond(1, nsites(lattice) + 1, 1, 1, 0), theta)
    @test_throws ArgumentError bond_phase(lattice, first(lattice.bonds), Inf)
    @test_throws ArgumentError gauge_angles(lattice, NaN)
end
