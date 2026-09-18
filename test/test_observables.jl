@testset "Cumulative transfer uses the measured profile" begin
    lattice = kagome_cylinder(3, 3)
    baseline = [0.15+0.002i for i in 1:27]
    current = copy(baseline)
    # More than one unit of cumulative spin must not be reduced modulo one.
    current[1:9] .-= 2/9
    current[19:27] .+= 2/9
    transfer = spin_transfer(lattice, baseline, current)
    @test length(transfer) == 2
    for t in transfer
        @test t.right ≈ 2 atol=1e-14
        @test t.left ≈ -2 atol=1e-14
        @test t.total ≈ 0 atol=1e-14
    end
    @test only(spin_transfer(lattice, current, baseline; cuts=[1])).right ≈ -2 atol=1e-14
    @test only(spin_transfer(lattice, baseline, baseline; cuts=[2])).right == 0
    @test_throws ArgumentError spin_transfer(lattice, baseline[1:26], current)
    @test_throws ArgumentError spin_transfer(lattice, fill(NaN, 27), current)
    @test_throws ArgumentError spin_transfer(lattice, baseline, current; cuts=[0])
    @test isempty(spin_transfer(kagome_cylinder(1, 3), zeros(9), zeros(9)))
end
