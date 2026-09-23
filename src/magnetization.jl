const _MAGNETIZATION_PRECISION_STATUSES = (
    "not_assessed", "not_run", "not_evaluated", "incomplete", "unmet", "passed")

_magnetization_require(condition, message) = condition || throw(ArgumentError(message))

function _magnetization_float(value, name)
    _magnetization_require(value isa Real && !(value isa Bool), "$name must be real, not Bool")
    number = try
        Float64(value)
    catch err
        err isa InterruptException && rethrow()
        throw(ArgumentError("$name must be representable as Float64"))
    end
    _magnetization_require(isfinite(number), "$name must be finite as Float64")
    return number
end

function _magnetization_charge(N, Q)
    _magnetization_require(Q isa Integer && !(Q isa Bool), "Q must be an integer, not Bool")
    return target_sector(N; Q).Q
end

# Arithmetic is exact for the normalized Float64 inputs, not for the unknown
# physical energies. Avoid cancellation in line comparisons and avoid using a
# tie tolerance to remove narrow but positive intervals.
_magnetization_exact(x::Float64) = Rational{BigInt}(x)

function _magnetization_minima(charges, energies, h, tolerance)
    values = [e - h * q / 2 for (q,e) in zip(charges,energies)]
    minimum_energy = minimum(values)
    minimizers = [q for (q,e) in zip(charges,values) if e == minimum_energy]
    nearby = [q for (q,e) in zip(charges,values) if e-minimum_energy <= tolerance]
    return (; minimum_energy, minimizers, nearby)
end

"""
    magnetization_curve(energies; N, fields=[], tie_atol=1e-10, precision_status=nothing)

Construct the zero-temperature lower envelope of `E0(Q) - h*Q/2` from total
zero-field sector energies. `energies` is a dictionary or an iterable of `Q=>E0`
pairs, with distinct integer `Q=2Sz` in `-N:2:N`. `M=Q/2`, `M/Msat=Q/N`.
Energy/site is not accepted as a different convention: callers must supply total
energies for the same model and geometry, with the uniform field removed.

Return a TOML-compatible schema-1 dictionary with sorted `sectors`, `coverage`,
all-sector `intervals`, lower-envelope `transitions`, and requested `samples`.
Interval status is `stable` (positive open width), `point` (one field only), or
`skipped` (never a minimizer). Missing competitors give unbounded endpoints;
coverage lists every missing Q and never assumes spin reversal. Completeness
means only that all charges were supplied, not that the energies are converged.

All real inputs are normalized to finite Float64. Exact rational comparisons
of these represented values determine topology and strict `minimizing_Q`.
`tie_atol` is an absolute total-energy tolerance used ONLY for sample
`near_minimizing_Q`. A sample is `tie` for multiple strict minimizers,
`near_tie` for additional nearby sectors, otherwise `unique`. It is not an
energy-error bound. Intervals and transitions never depend on this tolerance.

Finite bounds and sampled minimum energies are returned as Float64; overflow
or indistinguishable distinct bounds is rejected. Transition `field_is_exact`
states whether the exact crossing is representable as Float64. At a rounded
transition field, a sample may have a unique strict minimum despite the
transition's `coexisting_Q` at the exact crossing. Infinity denotes only a
genuinely unbounded interval endpoint.

Optional `precision_status` maps exactly the supplied Q values to strings
`not_assessed`, `not_run`, `not_evaluated`, `incomplete`, `unmet`, or `passed`.
These labels are copied, never inferred or certified. The function does not
load checkpoints, validate model/provenance, run DMRG, or establish a bulk
plateau. It reports the envelope of the supplied energies, including trials.
"""
function magnetization_curve(energies; N, fields=Float64[], tie_atol=1e-10,
        precision_status=nothing)
    _magnetization_require(N isa Integer && !(N isa Bool) && 0 < N < typemax(Int),
        "N must be a positive integer with N+1 representable as Int")
    N = Int(N)
    tolerance = _magnetization_float(tie_atol, "tie_atol")
    _magnetization_require(tolerance >= 0, "tie_atol must be nonnegative")
    _magnetization_require(fields isa AbstractVector || fields isa Tuple,
        "fields must be a vector, range, or tuple")
    requested_fields = [_magnetization_float(h, "field") for h in fields]
    _magnetization_require(applicable(iterate, energies), "energies must contain Q=>E0 pairs")
    table = Dict{Int,Float64}()
    for entry in energies
        _magnetization_require(entry isa Pair, "energies must contain Q=>E0 pairs")
        # Validate the scalar before sector arithmetic. In Julia 1.13, an
        # inlined target_sector before an always-throwing Bool/complex-energy
        # specialization can trap instead of propagating the ArgumentError.
        energy = _magnetization_float(last(entry), "sector energy")
        Q = _magnetization_charge(N, first(entry))
        _magnetization_require(!haskey(table,Q), "duplicate sector Q")
        table[Q] = energy
    end
    _magnetization_require(!isempty(table), "at least one sector energy is required")
    charges = sort!(collect(keys(table)))
    status = Dict(q => "not_assessed" for q in charges)
    if precision_status !== nothing
        _magnetization_require(precision_status isa AbstractDict,
            "precision_status must be a dictionary indexed by Q")
        status_charges = [_magnetization_charge(N,q) for q in keys(precision_status)]
        _magnetization_require(Set(status_charges) == Set(charges),
            "precision_status must describe exactly the supplied sectors")
        for (q,value) in precision_status
            _magnetization_require(value isa AbstractString && value in _MAGNETIZATION_PRECISION_STATUSES,
                "unknown sector precision status")
            status[q] = String(value)
        end
    end
    sectors = [Dict{String,Any}("Q" => q, "energy" => table[q], "M" => q/2,
        "magnetization" => q/N, "precision_status" => status[q]) for q in charges]
    exact_energies = [_magnetization_exact(table[q]) for q in charges]
    exact_tolerance = _magnetization_exact(tolerance)
    intervals = Dict{String,Any}[]
    breakpoints = Set{Rational{BigInt}}()
    for (i,q) in enumerate(charges)
        lower, upper = nothing, nothing
        for (j,other) in enumerate(charges)
            i == j && continue
            crossing = 2 * (exact_energies[j]-exact_energies[i]) / (big(other)-q)
            if other < q
                lower = lower === nothing ? crossing : max(lower,crossing)
            else
                upper = upper === nothing ? crossing : min(upper,crossing)
            end
        end
        label = lower === nothing || upper === nothing || lower < upper ? "stable" :
                lower == upper ? "point" : "skipped"
        lo = lower === nothing ? -Inf : _magnetization_float(lower,"lower bound")
        hi = upper === nothing ? Inf : _magnetization_float(upper,"upper bound")
        _magnetization_require(lower === nothing || upper === nothing ||
            lower == upper || lo != hi, "distinct interval bounds cannot be resolved as Float64")
        push!(intervals, Dict{String,Any}("Q" => q, "M" => q/2, "magnetization" => q/N,
            "h_lower" => lo, "h_upper" => hi, "status" => label))
        if label != "skipped"
            lower === nothing || push!(breakpoints,lower)
            upper === nothing || push!(breakpoints,upper)
        end
    end
    transitions = Dict{String,Any}[]
    for h in sort!(collect(breakpoints))
        exact = _magnetization_minima(charges,exact_energies,h,zero(exact_tolerance))
        field = _magnetization_float(h,"transition field")
        _magnetization_require(isempty(transitions) || last(transitions)["h"] < field,
            "distinct transitions cannot be resolved as Float64")
        push!(transitions, Dict{String,Any}("h" => field,
            "field_is_exact" => _magnetization_exact(field) == h,
            "from_Q" => first(exact.minimizers), "to_Q" => last(exact.minimizers),
            "coexisting_Q" => exact.minimizers,
            "coexisting_magnetization" => exact.minimizers ./ N))
    end
    samples = Dict{String,Any}[]
    for h in requested_fields
        values = _magnetization_minima(charges,exact_energies,_magnetization_exact(h),exact_tolerance)
        label = length(values.minimizers) > 1 ? "tie" : length(values.nearby) > 1 ? "near_tie" : "unique"
        push!(samples, Dict{String,Any}("h" => h,
            "minimum_energy" => _magnetization_float(values.minimum_energy,"minimum energy"),
            "minimizing_Q" => values.minimizers, "near_minimizing_Q" => values.nearby,
            "magnetization" => values.minimizers ./ N,
            "near_magnetization" => values.nearby ./ N, "status" => label))
    end
    missing = [q for q in -N:2:N if !haskey(table,q)]
    return Dict{String,Any}("schema_version" => 1, "format" => "KagomeDMRG.magnetization_curve",
        "N" => N, "scope" => "provided_sector_energies", "charge_convention" => "Q=2Sz",
        "energy_convention" => "total_zero_field_energy", "field_convention" => "E0-h*Q/2",
        "numerical_quality" => "not_assessed", "tie_atol" => tolerance,
        "precision_status_source" => precision_status === nothing ? "not_provided" : "caller",
        "sectors" => sectors, "coverage" => Dict("compared_Q" => charges,
            "missing_Q" => missing, "complete" => isempty(missing), "spin_reversal_assumed" => false),
        "intervals" => intervals, "transitions" => transitions, "samples" => samples)
end
