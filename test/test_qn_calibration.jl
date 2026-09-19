# Dense spin-basis calibration of QN truncation. The reference uses only
# LinearAlgebra on explicitly constructed coefficient matrices, never the
# production block decomposition, truncation helper, MPS, or MPO builders.
# These two/four-spin Q=0 controls are not kagome 1/9 plateau calculations.

function _qn_calibration_tensor(array, indices)
    shaped = reshape(array, Tuple(dim.(indices)))
    tensor = ITensor(ComplexF64, indices)
    for coordinate in CartesianIndices(shaped)
        value = shaped[coordinate]
        iszero(value) && continue
        tensor[(indices[j] => coordinate[j] for j in eachindex(indices))...] = value
    end
    return tensor
end

function _qn_calibration_fixture(probabilities; scale=1.0)
    length(probabilities) in (2, 4) || throw(ArgumentError("use two or four weights"))
    all(>=(0), probabilities) && sum(probabilities) > 0 ||
        throw(ArgumentError("weights must be nonnegative and nonzero"))
    n = length(probabilities) == 2 ? 2 : 4
    sites = [ITensors.siteind("S=1/2", j; conserve_qns=true) for j in 1:n]
    half = n ÷ 2
    d = 2^half
    matrix = zeros(ComplexF64, d, d)
    if n == 2
        matrix[1, 2] = sqrt(probabilities[1]) * cis(0.31)
        matrix[2, 1] = sqrt(probabilities[2]) * cis(-0.47)
    else
        # In Julia's column-major spin basis the two q=0 configurations are
        # rows/columns 2 and 3. Mix them with different complex unitaries on
        # the two sides so conjugation/transpose mistakes are observable.
        rotation(a, phase) = [cos(a) -cis(phase)*sin(a);
                              cis(-phase)*sin(a) cos(a)]
        ul = rotation(0.37, 0.63)
        ur = rotation(-0.51, -0.29)
        matrix[2:3, 2:3] = ul * Diagonal(sqrt.(probabilities[[1, 4]])) * ur'
        matrix[1, 4] = sqrt(probabilities[2]) * cis(0.79)
        matrix[4, 1] = sqrt(probabilities[3]) * cis(-0.43)
    end
    matrix .*= scale
    tensor = _qn_calibration_tensor(matrix, sites)
    charges = [half - 2count_ones(j - 1) for j in 1:d]
    return (; tensor, matrix, sites, left=sites[1:half], right=sites[(half + 1):end],
        charges, probabilities=Float64.(probabilities), scale=Float64(scale))
end

function _qn_calibration_density(fixture, side)
    side in (:left, :right) || throw(ArgumentError("unknown side"))
    m = fixture.matrix
    # The right reduced density has ket index first: transpose(M)*conj(M).
    # M'*M is its transpose and is different for these complex fixtures.
    density = side === :left ? m * m' : transpose(m) * conj(m)
    indices = side === :left ? fixture.left : fixture.right
    return (; density=Matrix(density), indices)
end

function _qn_calibration_dense_reference(density; maxdim, cutoff=0.0,
        use_absolute_cutoff=false, use_relative_cutoff=true)
    eigenvalues = reverse(LinearAlgebra.eigvals(Hermitian(density)))
    # All fixtures are positive semidefinite. Remove only roundoff on their
    # known null spaces, not a physical small-weight truncation threshold.
    floor = 100eps(Float64) * maximum(abs, eigenvalues)
    weights = [abs(value) <= floor ? 0.0 : max(value, 0.0) for value in eigenvalues]
    total = sum(weights)
    if use_absolute_cutoff
        required = max(1, count(>(cutoff), weights))
    else
        budget = cutoff * (use_relative_cutoff ? total : 1.0)
        required = something(findfirst(r -> sum(weights[(r + 1):end]) <= budget,
                                       1:length(weights)), length(weights))
    end
    rank = min(maxdim, required, length(weights))
    omitted = sum(weights[(rank + 1):end])
    error_scale = !use_absolute_cutoff && use_relative_cutoff ? total : 1.0
    return (; rank, weights, retained_weights=weights[1:rank], total,
        discarded_weight=omitted, discarded_probability=omitted / total,
        reported_error=omitted / error_scale, error_scale)
end

function _qn_calibration_metrics(fixture, reconstructed, isometry, link, spec,
        density, side; maxdim, cutoff=0.0,
        use_absolute_cutoff=false, use_relative_cutoff=true)
    indices = side === :left ? fixture.left : fixture.right
    d = size(fixture.matrix, 1)
    u = reshape(Array(isometry, indices..., link), d, dim(link))
    retained = reshape(Array(reconstructed, fixture.sites...), d, d)
    projector = u * u'
    projected = side === :left ? projector * fixture.matrix :
                                fixture.matrix * transpose(projector)
    reference = _qn_calibration_dense_reference(density; maxdim, cutoff,
        use_absolute_cutoff, use_relative_cutoff)
    input_norm2 = sum(abs2, fixture.matrix)
    density_norm = norm(density)
    selected = u' * density * u
    retained_charges = [dot(fixture.charges, abs2.(u[:, j])) for j in axes(u, 2)]
    charge_variances = [begin
        weights = abs2.(u[:, j])
        charge = retained_charges[j]
        dot((fixture.charges .- charge).^2, weights)
    end for j in axes(u, 2)]
    sector_error = sum(abs2(retained[i, j]) for i in 1:d, j in 1:d
                       if fixture.charges[i] + fixture.charges[j] != 0)
    return (; rank=dim(link), reference, reported_weights=collect(eigs(spec)),
        reported_error=truncerror(spec),
        retained_weights=sort(real.(eigvals(Hermitian(selected))); rev=true),
        wavefunction_loss=sum(abs2, fixture.matrix - retained) / input_norm2,
        projection_residual=norm(retained - projected) / sqrt(input_norm2),
        invariant_subspace_residual=norm(density * u - u * selected) / density_norm,
        orthogonality_error=norm(u' * u - I), retained_charges, charge_variances, sector_error,
        input_norm2)
end

"""Calibrate the public SVD directly, including either order of the partition."""
function _qn_calibrate_svd(fixture; side=:left, svd_function=svd, kwargs...)
    state = _qn_calibration_density(fixture, side)
    result = svd_function(fixture.tensor, state.indices; kwargs...)
    return _qn_calibration_metrics(fixture, result.U * result.S * result.V,
        result.U, result.u, result.spec, state.density, side; kwargs...)
end

"""Calibrate the public Hermitian QN eigen decomposition directly."""
function _qn_calibrate_eigen(fixture; side=:left, eigen_function=eigen, kwargs...)
    state = _qn_calibration_density(fixture, side)
    rows, columns = prime.(state.indices), dag.(state.indices)
    operator = _qn_calibration_tensor(state.density, [rows; columns])
    result = eigen_function(operator, rows, columns; ishermitian=true, kwargs...)
    isometry = replaceinds(result.Vt, rows, state.indices)
    # Project the original state without using the reported eigenvalues.
    projector = isometry * dag(prime(isometry, state.indices))
    projected = noprime(projector * prime(fixture.tensor, state.indices))
    return _qn_calibration_metrics(fixture, projected, isometry, result.l,
        result.spec, state.density, side; kwargs...)
end

"""
Calibrate factorize with an optional explicit PSD density-matrix perturbation.
For noise>0, the spectrum/error belong to rho+delta-rho; wavefunction_loss
still measures the projection loss of the original, unperturbed state.
"""
function _qn_calibrate_factorize(fixture; ortho="left", decomposition="svd",
        noise=0.0, factorize_function=factorize, kwargs...)
    side = ortho == "left" ? :left : :right
    state = _qn_calibration_density(fixture, side)
    density = copy(state.density)
    perturbation = nothing
    if noise > 0
        decomposition == "eigen" || throw(ArgumentError("noise requires eigen"))
        length(fixture.sites) == 4 || throw(ArgumentError("noise fixture needs four spins"))
        # This complex vector is confined to q=0 but does not commute with
        # that block's unperturbed density matrix. Its scale follows |psi|^2.
        direction = ComplexF64[0, 1, im, 0] / sqrt(2)
        delta = noise * sum(abs2, fixture.matrix) * (direction * direction')
        density += delta
        perturbation = _qn_calibration_tensor(delta,
            [prime.(state.indices); dag.(state.indices)])
    end
    left, right, spec, link = factorize_function(fixture.tensor, fixture.left;
        ortho, which_decomp=decomposition, eigen_perturbation=perturbation, kwargs...)
    isometry = side === :left ? left : right
    return _qn_calibration_metrics(fixture, left * right, isometry, link, spec,
        density, side; kwargs...)
end

function _qn_test_metrics(metrics; noiseless=true)
    ref = metrics.reference
    weight_tolerance = 3e-12 * ref.total
    @test metrics.rank == ref.rank
    @test length(metrics.reported_weights) == metrics.rank
    if length(metrics.reported_weights) == ref.rank
        @test sort(metrics.reported_weights; rev=true) ≈ ref.retained_weights atol=weight_tolerance rtol=0
    end
    if length(metrics.retained_weights) == ref.rank
        @test metrics.retained_weights ≈ ref.retained_weights atol=weight_tolerance rtol=0
    end
    @test metrics.reported_error ≈ ref.reported_error atol=3e-12*ref.total/ref.error_scale rtol=0
    @test metrics.projection_residual < 3e-12
    @test metrics.invariant_subspace_residual < 3e-12
    @test metrics.orthogonality_error < 3e-12
    @test all(<(3e-12), metrics.charge_variances)
    @test metrics.sector_error <= 1e-26 * metrics.input_norm2
    if noiseless
        @test metrics.wavefunction_loss ≈ ref.discarded_probability atol=3e-12 rtol=0
    end
    return nothing
end

function run_qn_calibration_tests(; svd_function=svd, eigen_function=eigen,
        factorize_function=factorize)
    @testset "Independent dense calibration of complex QN truncation" begin
        exact = [0.6, 0.2, 0.2, 0.0]
        near = [0.6, 0.2001, 0.1999, 0.0]
        # Explicit cases isolate the two repaired kernels, both complex
        # partitions, a three-way tie, and reconstruction through factorize.
        hard_cases = (
            (; name="exact SVD left", probabilities=exact, api=:svd, side=:left),
            (; name="exact eigen right", probabilities=exact, api=:eigen, side=:right),
            (; name="near SVD right", probabilities=near, api=:svd, side=:right),
            (; name="near eigen left", probabilities=near, api=:eigen, side=:left),
            (; name="three-way tie", probabilities=[0.3, 0.3, 0.1, 0.3], api=:svd, side=:left),
            (; name="nonzero tail", probabilities=[0.49, 0.24, 0.24, 0.03], api=:eigen, side=:right),
            (; name="factorize exact SVD", probabilities=exact, api=:factorize_svd, side=:right),
            (; name="factorize near eigen", probabilities=near, api=:factorize_eigen, side=:left))
        @testset "$(case.name)" for case in hard_cases
            fixture = _qn_calibration_fixture(case.probabilities)
            @test sort(svdvals(fixture.matrix).^2; rev=true) ≈
                  sort(case.probabilities; rev=true) atol=2e-13 rtol=0
            metrics = if case.api === :svd
                _qn_calibrate_svd(fixture; side=case.side, svd_function, maxdim=2, cutoff=0.0)
            elseif case.api === :eigen
                _qn_calibrate_eigen(fixture; side=case.side, eigen_function, maxdim=2, cutoff=0.0)
            else
                _qn_calibrate_factorize(fixture; ortho=String(case.side),
                    decomposition=case.api === :factorize_svd ? "svd" : "eigen",
                    factorize_function, maxdim=2, cutoff=0.0)
            end
            _qn_test_metrics(metrics)
        end

        # Each cutoff meaning is checked through both kernels at nonunit
        # norm. The relative rank-one case also distinguishes cutoff budgets.
        cutoff_cases = (
            (; name="relative sum SVD", api=:svd, scale=0.125,
               settings=(; maxdim=4, cutoff=0.205)),
            (; name="relative rank-one eigen", api=:eigen, scale=7.0,
               settings=(; maxdim=4, cutoff=0.405)),
            (; name="absolute sum SVD", api=:svd, scale=7.0,
               settings=(; maxdim=4, cutoff=0.205*7.0^2, use_relative_cutoff=false)),
            (; name="absolute sum eigen", api=:eigen, scale=0.125,
               settings=(; maxdim=4, cutoff=0.205*0.125^2, use_relative_cutoff=false)),
            (; name="absolute value SVD", api=:svd, scale=0.125,
               settings=(; maxdim=4, cutoff=0.205*0.125^2, use_absolute_cutoff=true)),
            (; name="absolute value eigen", api=:eigen, scale=7.0,
               settings=(; maxdim=4, cutoff=0.205*7.0^2, use_absolute_cutoff=true)))
        @testset "$(case.name)" for case in cutoff_cases
            fixture = _qn_calibration_fixture(near; scale=case.scale)
            metrics = case.api === :svd ?
                _qn_calibrate_svd(fixture; svd_function, case.settings...) :
                _qn_calibrate_eigen(fixture; eigen_function, case.settings...)
            _qn_test_metrics(metrics)
        end

        noise_cases = ((; ortho="left", scale=1.0), (; ortho="right", scale=7.0))
        @testset "Perturbed density ortho=$(case.ortho) scale=$(case.scale)" for case in noise_cases
            fixture = _qn_calibration_fixture(exact; scale=case.scale)
            metrics = _qn_calibrate_factorize(fixture; ortho=case.ortho, decomposition="eigen",
                noise=0.1, factorize_function, maxdim=2, cutoff=0.0)
            _qn_test_metrics(metrics; noiseless=false)
            @test abs(metrics.wavefunction_loss - metrics.reported_error) > 1e-3
        end

        @testset "Per-sector minima share the hard global budget" begin
            # The mandatory q=0,+2,-2 states have weights 0.6,0.3,0.09.
            # Their analytic constrained optimum needs three slots even when
            # cutoff=1 would otherwise keep one. The remaining loss is 0.01.
            fixture = _qn_calibration_fixture([0.6, 0.3, 0.09, 0.01])
            state = _qn_calibration_density(fixture, :left)
            result = svd_function(fixture.tensor, state.indices;
                min_blockdim=1, maxdim=3, cutoff=1.0)
            metrics = _qn_calibration_metrics(fixture,
                result.U * result.S * result.V, result.U, result.u,
                result.spec, state.density, :left; maxdim=3, cutoff=0.0)
            _qn_test_metrics(metrics)
            @test sort(metrics.retained_charges) ≈ [-2.0, 0.0, 2.0] atol=3e-12 rtol=0
            @test_throws ArgumentError svd_function(fixture.tensor, state.indices;
                min_blockdim=1, maxdim=2, cutoff=0.0)
        end
    end
    return nothing
end

run_qn_calibration_tests()
