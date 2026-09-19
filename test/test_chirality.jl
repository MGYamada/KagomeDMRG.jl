# Independent references use Cartesian Pauli components and periodic Euclidean
# distances, rather than the production ladder expansion or triangle template.
function _chirality_pauli_reference()
    spin = (ComplexF64[0 1; 1 0] / 2,
            ComplexF64[0 -im; im 0] / 2,
            ComplexF64[1 0; 0 -1] / 2)
    result = zeros(ComplexF64, 8, 8)
    for a in 1:3, b in 1:3, c in 1:3
        epsilon = (a - b) * (b - c) * (c - a) / 2
        # Array(ITensor, sites...) makes physical site 1 the fastest index.
        result .+= epsilon .* kron(spin[c], spin[b], spin[a])
    end
    return result
end

function _chirality_triangle_key(triangle)
    order = sortperm(collect(triangle.sites))
    anchor = triangle.image_y[first(order)]
    return (ntuple(a -> triangle.sites[order[a]], 3),
            ntuple(a -> triangle.image_y[order[a]] - anchor, 3))
end

function _chirality_geometric_triangles(lattice)
    expected = Set{Tuple{NTuple{3,Int},NTuple{3,Int}}}()
    wrap = (lattice.Ly / 2, sqrt(3.0) * lattice.Ly / 2)
    lifted(i, m) = ntuple(d -> lattice.sites[i].position[d] + m * wrap[d], 2)
    near(a, b) = isapprox(sum((a[d] - b[d])^2 for d in 1:2), 1 / 4;
                          atol=2e-14, rtol=0)
    # Anchor the smallest index at image zero. A third edge must close in the
    # same plane, excluding a spurious loop that winds around the cylinder.
    N = nsites(lattice)
    for i in 1:(N - 2), j in (i + 1):(N - 1), mj in -1:1
        ri, rj = lifted(i, 0), lifted(j, mj)
        near(ri, rj) || continue
        for k in (j + 1):N, mk in -1:1
            rk = lifted(k, mk)
            near(ri, rk) && near(rj, rk) &&
                push!(expected, ((i, j, k), (0, mj, mk)))
        end
    end
    return expected
end

function _chirality_three_spin_mps(sites, amplitudes)
    tensor = ITensor(ComplexF64, sites...)
    for down in 1:3
        tensor[(sites[a] => (a == down ? 2 : 1) for a in 1:3)...] = amplitudes[down]
    end
    return MPS(tensor, sites; cutoff=0.0)
end

function _chirality_embedded_mps(sites, vertices, amplitudes, angles)
    # Exactly three product components, Q=2 on N=18, with no dense 18-site
    # tensor. This synthetic chiral control is not a ground-state claim.
    @assert length(sites) == 18
    spectators = setdiff(1:18, collect(vertices))
    pieces = MPS[]
    for down in 1:3
        labels = fill("Dn", 18)
        labels[spectators[1:8]] .= "Up"
        labels[collect(vertices)] .= "Up"
        labels[vertices[down]] = "Dn"
        phase = sum(angles[i] * (labels[i] == "Up" ? 0.5 : -0.5) for i in 1:18)
        push!(pieces, amplitudes[down] * cis(phase) * MPS(ComplexF64, sites, labels))
    end
    return add(add(pieces[1], pieces[2]; alg="directsum"), pieces[3]; alg="directsum")
end

@testset "Oriented scalar chirality calibration" begin
    metrics = Dict{Symbol,Any}()
    @testset "Independent three-spin Cartesian reference" begin
        sites = siteinds("S=1/2", 3; conserve_qns=true)
        reference = _chirality_pauli_reference()
        C = scalar_chirality_mpo(sites, (1, 2, 3))
        full = full_mpo_matrix(C, sites)
        metrics[:matrix_max_error] = maximum(abs.(full - reference))
        @test full ≈ reference atol=2e-14 rtol=0
        @test full ≈ full' atol=2e-14 rtol=0
        charge = Diagonal([3 - 2count_ones(bits) for bits in 0:7])
        @test norm(charge * full - full * charge) < 2e-14
        spectrum = eigvals(Hermitian(full))
        metrics[:eigenvalues] = spectrum
        @test spectrum ≈ [-sqrt(3) / 4, -sqrt(3) / 4, 0, 0, 0, 0,
                           sqrt(3) / 4, sqrt(3) / 4] atol=2e-14 rtol=0
        @test full_mpo_matrix(scalar_chirality_mpo(sites, (2, 3, 1)), sites) ≈
              full atol=2e-14 rtol=0
        @test full_mpo_matrix(scalar_chirality_mpo(sites, (1, 3, 2)), sites) ≈
              -full atol=2e-14 rtol=0

        omega = cis(2pi / 3)
        amplitudes = ComplexF64[1, omega, omega^2] / sqrt(3)
        vector = zeros(ComplexF64, 8)
        vector[[2, 3, 5]] = amplitudes
        @test reference * vector ≈ sqrt(3) / 4 * vector atol=2e-14 rtol=0
        @test reference * conj(vector) ≈ -sqrt(3) / 4 * conj(vector) atol=2e-14 rtol=0
        psi = _chirality_three_spin_mps(sites, amplitudes)
        known_value = real(inner(psi', C, psi))
        metrics[:known_positive_chirality] = known_value
        @test flux(psi) == QN("Sz", 1)
        @test known_value ≈ sqrt(3) / 4 atol=2e-14 rtol=0
        product = MPS(ComplexF64, sites, fill("Up", 3))
        @test abs(inner(product', C, product)) < 2e-14

        angles = (0.27, -0.41, 0.63)
        rotations = Diagonal([cis(sum(angles[a] * (isodd(bits >> (a - 1)) ? -0.5 : 0.5)
                                       for a in 1:3)) for bits in 0:7])
        rotated = full_mpo_matrix(scalar_chirality_mpo(sites, (1, 2, 3); angles), sites)
        @test rotated ≈ rotations * reference * rotations' atol=2e-14 rtol=0
        @test rotated ≈ rotated' atol=2e-14 rtol=0
        @test_throws ArgumentError scalar_chirality_mpo(sites, (1, 2, 4))
        @test_throws ArgumentError scalar_chirality_mpo(sites, (1, 2, 3); angles=(0, NaN, 0))
        @test_throws ArgumentError scalar_chirality_mpo(siteinds("S=1/2", 3), (1, 2, 3))
    end

    @testset "Independent periodic triangles and CCW orientation" begin
        # One-column termination and a distinct circumference with down/seam
        # triangles. The N18 measurement below independently exercises Ly=3.
        for (lx, ly) in ((1, 3), (2, 4))
            lattice = kagome_cylinder(lx, ly)
            triangles = oriented_triangles(lattice)
            keys = _chirality_triangle_key.(triangles)
            @test length(triangles) == (2lx - 1) * ly
            @test length(Set(keys)) == length(keys)
            @test Set(keys) == _chirality_geometric_triangles(lattice)
            @test count(t -> t.kind === :up, triangles) == lx * ly
            @test count(t -> t.kind === :down, triangles) == (lx - 1) * ly
            wrap = (ly / 2, sqrt(3.0) * ly / 2)
            for triangle in triangles
                p = ntuple(a -> ntuple(d -> lattice.sites[triangle.sites[a]].position[d] +
                                      triangle.image_y[a] * wrap[d], 2), 3)
                area2 = (p[2][1] - p[1][1]) * (p[3][2] - p[1][2]) -
                        (p[2][2] - p[1][2]) * (p[3][1] - p[1][1])
                @test area2 ≈ sqrt(3) / 8 atol=2e-14 rtol=0
                bottom = minimum(point[2] for point in p)
                geometric_kind = count(point -> isapprox(point[2], bottom; atol=2e-14), p) == 2 ?
                                 :up : :down
                @test triangle.kind === geometric_kind
                for (a, b) in ((1, 2), (2, 3), (3, 1))
                    i, j = triangle.sites[a], triangle.sites[b]
                    bond = only(filter(edge -> (edge.i == i && edge.j == j) ||
                                                (edge.i == j && edge.j == i), lattice.bonds))
                    winding = bond.i == i ? bond.wy : -bond.wy
                    @test triangle.image_y[b] - triangle.image_y[a] == winding
                end
            end
        end
        # Elementary NN plaquettes are geometric observables even when the
        # Hamiltonian contains additional J2/J3 bonds and changed couplings.
        nn = oriented_triangles(kagome_cylinder(2, 3))
        extended = oriented_triangles(kagome_j1j2j3_cylinder(2, 3; J1=0.7, J2=0.5, J3=-0.2))
        triangle_data(t) = (t.sites, t.image_y, t.kind)
        @test triangle_data.(extended) == triangle_data.(nn)
        @test_throws ArgumentError KagomeTriangle((1, 1, 3), (0, 0, 0), :up)
        @test_throws ArgumentError KagomeTriangle((1, 2, 3), (0, 0, 0), :invalid)
    end

    @testset "N18 fixed-charge seam control and gauge covariance" begin
        lattice = kagome_cylinder(2, 3)
        sites = spin_sites(lattice)
        triangles = oriented_triangles(lattice)
        # Independent seam-crossing down triangle: A(1,0), B(0,0), C(1,2).
        vertices = (10, 2, 18)
        selected = only(findall(t -> Set(t.sites) == Set(vertices), triangles))
        theta = 0.61
        seam_angles = zeros(18)
        seam_angles[18] = theta  # C image -1 gives alpha_C = +theta.
        uniform_angles = [-theta * (s.y + (s.sublattice === :C ? 0.5 : 0.0)) / lattice.Ly
                          for s in lattice.sites]
        omega = cis(2pi / 3)
        amplitudes = ComplexF64[1, omega, omega^2] / sqrt(3)
        psi = _chirality_embedded_mps(sites, vertices, amplitudes, seam_angles)
        rotated = _chirality_embedded_mps(sites, vertices, amplitudes, seam_angles + uniform_angles)
        @test flux(psi) == flux(rotated) == QN("Sz", 2)
        @test maxlinkdim(psi) <= 3
        psi[1] *= 2cis(0.19)
        before = deepcopy(psi)
        values = triangle_chiralities(psi, lattice, theta)
        expected = zeros(length(triangles))
        expected[selected] = sqrt(3) / 4
        @test eltype(values) <: Real
        @test values ≈ expected atol=2e-12 rtol=0
        uniform = triangle_chiralities(rotated, lattice, theta; gauge=:uniform)
        gauge_error = maximum(abs.(uniform - values))
        preservation_error = maximum(norm(psi[i] - before[i]) for i in 1:18)
        metrics[:seam_uniform_max_error] = gauge_error
        metrics[:normalized_control_max_error] = maximum(abs.(values - expected))
        metrics[:input_tensor_preservation_error] = preservation_error
        @test uniform ≈ values atol=2e-12 rtol=0
        @test preservation_error == 0
        @test norm(psi) ≈ 2 atol=2e-13 rtol=0
        for (a, b) in ((1, 2), (2, 3), (3, 1))
            i, j = vertices[a], vertices[b]
            bond = only(filter(edge -> (edge.i == i && edge.j == j) ||
                                        (edge.i == j && edge.j == i), lattice.bonds))
            orientation = bond.i == i ? 1 : -1
            @test cis(seam_angles[i] - seam_angles[j]) ≈
                  cis(orientation * bond_phase(lattice, bond, theta)) atol=2e-14
            @test cis(seam_angles[i] + uniform_angles[i] - seam_angles[j] - uniform_angles[j]) ≈
                  cis(orientation * bond_phase(lattice, bond, theta; gauge=:uniform)) atol=2e-14
        end
        @test_throws ArgumentError triangle_chiralities(psi, lattice, NaN)
        @test_throws ArgumentError triangle_chiralities(psi, lattice, theta; gauge=:invalid)
        zero = deepcopy(psi)
        zero[1] *= 0
        @test_throws ArgumentError triangle_chiralities(zero, lattice)
    end
    @info "Independent chirality calibration" metrics
end
