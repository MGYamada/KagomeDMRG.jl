using TOML

include(joinpath(@__DIR__, "..", "examples", "magnetization_curve.jl"))

@testset "Finite-system magnetization envelope" begin
    # Exact hand calculation: E(Q)=Q^2/4 for every allowed charge of N=4.
    # Unsorted inputs and fields exercise deterministic sector order while
    # retaining the requested sample order (including repeated fields).
    energies = [4 => 4, -2 => 1, 0 => 0, -4 => 4, 2 => 1]
    fields = [4, -4, 0, 1, 1.25, 1]
    curve = magnetization_curve(energies; N=4, fields)
    @test curve["format"] == "KagomeDMRG.magnetization_curve"
    @test curve["schema_version"] == 1
    @test [row["Q"] for row in curve["sectors"]] == [-4, -2, 0, 2, 4]
    @test [row["energy"] for row in curve["sectors"]] == [4, 1, 0, 1, 4]
    @test [row["M"] for row in curve["sectors"]] == [-2, -1, 0, 1, 2]
    @test [row["magnetization"] for row in curve["sectors"]] == [-1, -0.5, 0, 0.5, 1]
    @test all(row["energy"] isa Float64 for row in curve["sectors"])
    @test all(row["precision_status"] == "not_assessed" for row in curve["sectors"])
    @test curve["coverage"]["compared_Q"] == [-4, -2, 0, 2, 4]
    @test isempty(curve["coverage"]["missing_Q"])
    @test curve["coverage"]["complete"]
    @test !curve["coverage"]["spin_reversal_assumed"]

    intervals = curve["intervals"]
    @test [row["Q"] for row in intervals] == [-4, -2, 0, 2, 4]
    @test [row["h_lower"] for row in intervals] == [-Inf, -3, -1, 1, 3]
    @test [row["h_upper"] for row in intervals] == [-3, -1, 1, 3, Inf]
    @test all(row["status"] == "stable" for row in intervals)
    @test [row["M"] for row in intervals] == [-2, -1, 0, 1, 2]
    @test [row["magnetization"] for row in intervals] == [-1, -0.5, 0, 0.5, 1]
    transitions = curve["transitions"]
    @test [row["h"] for row in transitions] == [-3, -1, 1, 3]
    @test [row["from_Q"] for row in transitions] == [-4, -2, 0, 2]
    @test [row["to_Q"] for row in transitions] == [-2, 0, 2, 4]
    @test [row["coexisting_Q"] for row in transitions] == [[-4, -2], [-2, 0], [0, 2], [2, 4]]
    samples = curve["samples"]
    @test [row["h"] for row in samples] == fields
    @test all(row["h"] isa Float64 for row in samples)
    @test [row["minimum_energy"] for row in samples] == [-4, -4, 0, 0, -0.25, 0]
    @test [row["minimizing_Q"] for row in samples] == [[4], [-4], [0], [0, 2], [2], [0, 2]]
    @test [row["near_minimizing_Q"] for row in samples] == [row["minimizing_Q"] for row in samples]
    @test [row["status"] for row in samples] == ["unique", "unique", "unique", "tie", "unique", "tie"]

    # Infinity denotes an unbounded interval and survives the public TOML record.
    buffer = IOBuffer()
    TOML.print(buffer, curve; sorted=true)
    @test TOML.parse(String(take!(buffer))) == curve
    @test magnetization_curve((q => e for (q, e) in energies); N=4, fields) == curve
end

@testset "Energy-table example and tie-preserving output" begin
    input = joinpath(@__DIR__, "..", "examples", "configs", "magnetization_curve_demo.toml")
    config = TOML.parsefile(input)
    @test_throws ArgumentError MagnetizationCurveExample.main(String[])
    @test_throws ArgumentError MagnetizationCurveExample.analyze_config(
        merge(config, Dict("unknown_option" => true)))
    labeled = deepcopy(config)
    labeled["sectors"][1]["precision_status"] = "unmet"
    labels = MagnetizationCurveExample.analyze_config(labeled)["curve"]["sectors"]
    @test first(labels)["precision_status"] == "unmet"
    @test all(row["precision_status"] == "not_assessed" for row in labels[2:end])
    mktempdir() do directory
        output = joinpath(directory, "demo")
        report = MagnetizationCurveExample.analyze(output, input)
        @test report["input_scope"] == "synthetic_fixture"
        @test TOML.parsefile(joinpath(output, "report.toml")) == report
        csv = read(joinpath(output, "samples.csv"), String)
        @test length(split(chomp(csv), '\n')) == 10
        @test occursin("1.0,0.0,\"tie\",\"0;2\",\"0;2\",\"0.0;0.5\",\"0.0;0.5\"", csv)
        before = _checkpoint_test_file_hashes(output)
        @test_throws ArgumentError MagnetizationCurveExample.analyze(output, input)
        @test _checkpoint_test_file_hashes(output) == before
        invalid = joinpath(directory, "invalid.toml")
        open(invalid, "w") do io
            TOML.print(io, merge(config, Dict("schema_version" => true)))
        end
        rejected = joinpath(directory, "rejected")
        @test_throws ArgumentError MagnetizationCurveExample.analyze(rejected, invalid)
        @test !ispath(rejected)
    end
end

@testset "Skipped sectors, point coexistence, and genuine narrow intervals" begin
    jump = magnetization_curve(Dict(-2 => 0, 0 => 1, 2 => 0); N=2, fields=[-1, 0, 1])
    intervals = Dict(row["Q"] => row for row in jump["intervals"])
    @test intervals[0]["status"] == "skipped"
    @test intervals[0]["h_lower"] == 1
    @test intervals[0]["h_upper"] == -1
    transition = only(jump["transitions"])
    @test transition["h"] == 0
    @test transition["from_Q"] == -2
    @test transition["to_Q"] == 2
    @test transition["coexisting_Q"] == [-2, 2]
    @test [row["minimizing_Q"] for row in jump["samples"]] == [[-2], [-2, 2], [2]]

    point = magnetization_curve([-2 => 0, 0 => 0, 2 => 0]; N=2, fields=[0])
    @test [row["status"] for row in point["intervals"]] == ["stable", "point", "stable"]
    @test point["intervals"][2]["h_lower"] == point["intervals"][2]["h_upper"] == 0
    @test only(point["transitions"])["coexisting_Q"] == [-2, 0, 2]
    @test only(point["samples"])["minimizing_Q"] == [-2, 0, 2]
    @test only(point["samples"])["status"] == "tie"

    # A tolerance for sampled energies must never merge two real transitions.
    delta = 2.0^-40
    narrow = magnetization_curve([-2 => delta, 0 => 0, 2 => delta];
        N=2, fields=[0], tie_atol=1e-6)
    middle = narrow["intervals"][2]
    @test middle["status"] == "stable"
    @test middle["h_lower"] == -delta
    @test middle["h_upper"] == delta
    @test [row["h"] for row in narrow["transitions"]] == [-delta, delta]
    sample = only(narrow["samples"])
    @test sample["minimizing_Q"] == [0]
    @test sample["near_minimizing_Q"] == [-2, 0, 2]
    @test sample["status"] == "near_tie"
    strict = magnetization_curve([-2 => delta, 0 => 0, 2 => delta];
        N=2, fields=[0], tie_atol=0)
    @test strict["intervals"] == narrow["intervals"]
    @test strict["transitions"] == narrow["transitions"]
    @test only(strict["samples"])["near_minimizing_Q"] == [0]
    @test only(strict["samples"])["status"] == "unique"

    # Absolute (not relative) tolerance is invariant under this exactly
    # representable common energy offset, even far from the origin.
    base = magnetization_curve([-2 => 0, 0 => 0.125, 2 => 1];
        N=2, fields=[0], tie_atol=0.25)
    shifted = magnetization_curve([-2 => 2.0^40, 0 => 2.0^40 + 0.125, 2 => 2.0^40 + 1];
        N=2, fields=[0], tie_atol=0.25)
    @test shifted["intervals"] == base["intervals"]
    @test shifted["transitions"] == base["transitions"]
    @test only(base["samples"])["near_minimizing_Q"] == [-2, 0]
    @test only(shifted["samples"])["near_minimizing_Q"] == [-2, 0]
    @test only(shifted["samples"])["minimizing_Q"] == [-2]
    @test only(shifted["samples"])["status"] == "near_tie"

    # A Float64 sample energy may round away differences which still choose
    # different sectors. Strict minimization uses represented inputs exactly.
    cancellation = magnetization_curve([-2 => 1e16, 0 => 1e16, 2 => 1e16];
        N=2, fields=[1], tie_atol=0.25)
    @test only(cancellation["samples"])["minimizing_Q"] == [2]
    @test only(cancellation["samples"])["near_minimizing_Q"] == [2]
    @test only(cancellation["samples"])["status"] == "unique"

    # 1/3 cannot be represented exactly as Float64. Coexistence belongs to the
    # exact transition; sampling its rounded-down display value is not a tie.
    rounded = magnetization_curve([-6 => 0, 0 => 1]; N=6,
        fields=[1/3], tie_atol=0)
    @test only(rounded["transitions"])["h"] == 1/3
    @test !only(rounded["transitions"])["field_is_exact"]
    @test only(rounded["transitions"])["coexisting_Q"] == [-6, 0]
    @test only(rounded["samples"])["minimizing_Q"] == [-6]
    @test only(rounded["samples"])["status"] == "unique"
end

@testset "Partial coverage and detached precision metadata" begin
    energies = Dict(0 => 0.0, 2 => 1.0)
    statuses = Dict(0 => "unmet", 2 => "passed")
    fields = [0.0, 2.0]
    curve = magnetization_curve(energies; N=4, fields, precision_status=statuses)
    @test curve["coverage"]["compared_Q"] == [0, 2]
    @test curve["coverage"]["missing_Q"] == [-4, -2, 4]
    @test !curve["coverage"]["complete"]
    @test !curve["coverage"]["spin_reversal_assumed"]
    @test [row["precision_status"] for row in curve["sectors"]] == ["unmet", "passed"]
    @test only(curve["transitions"])["h"] == 1
    @test [row["minimizing_Q"] for row in curve["samples"]] == [[0], [2]]
    # Mutating either side must not rewrite the other or synthesize symmetry partners.
    energies[0] = -100
    statuses[2] = "incomplete"
    fields[1] = 100
    @test curve["sectors"][1]["energy"] == 0
    @test curve["sectors"][2]["precision_status"] == "passed"
    @test curve["samples"][1]["h"] == 0
    curve["coverage"]["compared_Q"][1] = -4
    @test curve["sectors"][1]["Q"] == 0
    @test Set(keys(energies)) == Set([0, 2])

    singleton = magnetization_curve([1 => -2//3]; N=3, fields=[-100, 100])
    @test singleton["coverage"]["missing_Q"] == [-3, -1, 3]
    @test !singleton["coverage"]["complete"]
    @test isempty(singleton["transitions"])
    @test only(singleton["intervals"])["h_lower"] == -Inf
    @test only(singleton["intervals"])["h_upper"] == Inf
    @test only(singleton["intervals"])["status"] == "stable"
    @test only(singleton["sectors"])["energy"] === Float64(-2//3)
    @test all(row["minimizing_Q"] == [1] for row in singleton["samples"])
    @test isempty(magnetization_curve([1 => 0]; N=1)["samples"])

    labels = ["not_assessed", "not_run", "not_evaluated", "incomplete", "unmet", "passed"]
    metadata = Dict(q => status for (q, status) in zip(-5:2:5, labels))
    retained = magnetization_curve(Dict(q => 0 for q in -5:2:5);
        N=5, precision_status=metadata)
    @test [row["precision_status"] for row in retained["sectors"]] == labels
end

@testset "Magnetization analysis input validation" begin
    valid = [-2 => 0, 0 => 0, 2 => 0]
    for bad_N in (0, -1, true, 2.0)
        @test_throws ArgumentError magnetization_curve(valid; N=bad_N)
    end
    for bad_energies in (Pair{Int,Float64}[], [-2 => 0, -2 => 1], [0 => 0, 1 => 1],
                         [-4 => 0], [0.0 => 0], [false => 0], [0 => NaN], [0 => Inf],
                         [0 => true], [0 => 1 + 0im], [(0, 0)])
        @test_throws ArgumentError magnetization_curve(bad_energies; N=2)
    end
    for bad_fields in ([NaN], [Inf], [-Inf], [true], [1 + 0im])
        @test_throws ArgumentError magnetization_curve(valid; N=2, fields=bad_fields)
    end
    for bad_tolerance in (-1, NaN, Inf, true, 1 + 0im)
        @test_throws ArgumentError magnetization_curve(valid; N=2, tie_atol=bad_tolerance)
    end
    for bad_status in (Dict(0 => "passed"), Dict(-2 => "passed", 0 => "passed", 2 => "bad"),
                       Dict(-2 => "passed", 0 => "passed", 2 => true),
                       Dict(-2.0 => "passed", 0.0 => "passed", 2.0 => "passed"),
                       Dict(-2 => "passed", 0 => "passed", 2 => "passed", 4 => "passed"))
        @test_throws ArgumentError magnetization_curve(valid; N=2, precision_status=bad_status)
    end
    # Refuse to serialize two distinct exact bounds as the same Float64, or
    # overflow a finite crossing/minimum into a spurious infinite result.
    @test_throws ArgumentError magnetization_curve(
        [-4 => nextfloat(0.0), 0 => 0.0, 4 => nextfloat(0.0)]; N=4)
    @test_throws ArgumentError magnetization_curve(
        [0 => floatmax(Float64), 2 => -floatmax(Float64)]; N=2)
    @test_throws ArgumentError magnetization_curve(
        [2 => -floatmax(Float64)]; N=2, fields=[floatmax(Float64)])
end

@testset "Magnetization curve against a full spin-basis field Hamiltonian" begin
    # Bounded ED only: ten sectors of dimension <=126 plus three 512x512
    # full-space solves, with BLAS=1 set by runtests.jl. This reference builds
    # bonds by Cartesian distance and never calls a production MPO builder.
    N = 9
    energies = Dict(q => eigmin(Hermitian(reference_hamiltonian(1, 3; Q=q).H))
                    for q in -N:2:N)
    fields = [-0.2, 0.2, 4.0]
    curve = magnetization_curve(energies; N, fields)
    @test curve["coverage"]["complete"]
    @test energies[-N] ≈ 3 atol=1e-13 rtol=0
    @test energies[N] ≈ 3 atol=1e-13 rtol=0

    full = reference_hamiltonian(1, 3; nup=:all)
    # Independently sum physical site spins from the bit representation.
    total_sz = [sum(iszero(state & (UInt64(1) << (i - 1))) ? -0.5 : 0.5
                    for i in 1:N) for state in full.basis]
    for (h, sample) in zip(fields, curve["samples"])
        field_ground = eigen(Hermitian(full.H - h * Diagonal(total_sz)))
        psi = field_ground.vectors[:, 1]
        measured_M = real(dot(psi, total_sz .* psi))
        @test sample["status"] == "unique"
        @test sample["minimum_energy"] ≈ field_ground.values[1] atol=2e-12 rtol=0
        @test only(sample["minimizing_Q"])/2 ≈ measured_M atol=2e-12 rtol=0
        @test only(sample["minimizing_Q"])/N ≈ 2measured_M/N atol=5e-13 rtol=0
    end
end
