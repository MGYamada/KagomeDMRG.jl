# Select block ranks from the same globally ordered entries used by Spectrum.
# A scalar cutoff threshold cannot select exactly n states at a tied boundary.
# Equal weights are ordered by block coordinate, then by within-block index.
# Within each block weights must already be sorted in descending order, as
# provided by the dense SVD/eigen kernels. Each retained block is then a prefix.
function _truncate_block_spectrum(
        block_weights::AbstractVector{<:AbstractVector}, block_keys;
        min_blockdim = nothing,
        mindim = nothing,
        maxdim = nothing,
        cutoff = nothing,
        use_absolute_cutoff = nothing,
        use_relative_cutoff = nothing
    )
    length(block_weights) == length(block_keys) ||
        throw(ArgumentError("block weights and coordinates must have equal length"))
    isempty(block_weights) && throw(ArgumentError("cannot truncate an empty block spectrum"))
    minimum_block = replace_nothing(min_blockdim, 0)
    minimum_block isa Integer && minimum_block >= 0 ||
        throw(ArgumentError("min_blockdim must be a nonnegative integer"))

    weights = eltype(first(block_weights))[]
    block_ids = Int[]
    local_ids = Int[]
    # Stable global sorting now has a canonical order for exact ties, even if
    # the input stores nonzero blocks in a different insertion order.
    for b in sortperm(block_keys)
        append!(weights, block_weights[b])
        append!(block_ids, fill(b, length(block_weights[b])))
        append!(local_ids, eachindex(block_weights[b]))
    end
    isempty(weights) && throw(ArgumentError("cannot truncate an empty block spectrum"))
    order = sortperm(weights; rev = true, alg = Base.Sort.MergeSort)
    sorted_weights = weights[order]

    selected = falses(length(weights))
    if iszero(minimum_block)
        # Preserve the upstream cutoff, mindim, normalization and hard maxdim
        # rules. Only discard the threshold: the chosen rank is authoritative.
        truncerr, _ = truncate!(
            sorted_weights; mindim, maxdim, cutoff,
            use_absolute_cutoff, use_relative_cutoff
        )
        selected[order[1:length(sorted_weights)]] .= true
    else
        # Nonzero per-block minima are constraints, not extra states added
        # after truncation. Fill the remaining hard-maxdim slots by weight.
        capacity = min(replace_nothing(maxdim, length(weights)), length(weights))
        capacity isa Integer && capacity >= 1 ||
            throw(ArgumentError("maxdim must permit at least one retained state"))
        mandatory = local_ids .<= minimum_block
        selected .= mandatory
        nkeep = count(selected)
        nkeep <= capacity || throw(ArgumentError(
            "min_blockdim requires $nkeep states, exceeding maxdim=$capacity"
        ))
        for i in order
            nkeep == capacity && break
            selected[i] && continue
            selected[i] = true
            nkeep += 1
        end

        absolute = replace_nothing(use_absolute_cutoff, default_use_absolute_cutoff(weights))
        relative = replace_nothing(use_relative_cutoff, default_use_relative_cutoff(weights))
        threshold = replace_nothing(cutoff, typemin(eltype(weights)))
        minimum_total = max(1, replace_nothing(mindim, default_mindim(weights)))
        scale = !absolute && relative ? sum(sorted_weights) : one(eltype(weights))
        iszero(scale) && (scale = one(eltype(weights)))
        discarded = zero(eltype(weights))
        for i in Iterators.reverse(order)
            selected[i] || (discarded += weights[i])
        end
        # Discard the least important eligible state first. Once its removal
        # would exceed the cutoff, no larger eligible weight can be discarded.
        for i in Iterators.reverse(order)
            nkeep <= minimum_total && break
            (!selected[i] || mandatory[i]) && continue
            allowed = absolute ? weights[i] <= threshold :
                discarded + weights[i] <= threshold * scale
            allowed || break
            selected[i] = false
            nkeep -= 1
            discarded += weights[i]
        end
        sorted_weights = weights[filter(i -> selected[i], order)]
        truncerr = discarded / scale
    end

    blockdims = zeros(Int, length(block_weights))
    for i in eachindex(selected)
        selected[i] && (blockdims[block_ids[i]] += 1)
    end
    return blockdims, sorted_weights, truncerr
end
