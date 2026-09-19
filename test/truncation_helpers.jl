# Shared MPS integration controls. The independent dense calibration in
# test_qn_calibration.jl intentionally does not use this MPO/MPS fixture.
function _heisenberg_chain_control(n; staggered_field=0.0)
    sites = siteinds("S=1/2", n; conserve_qns=true)
    terms = OpSum()
    for j in 1:(n-1)
        terms += "Sz", j, "Sz", j+1
        terms += 0.5, "S+", j, "S-", j+1
        terms += 0.5, "S-", j, "S+", j+1
    end
    if !iszero(staggered_field)
        for j in 1:n
            terms += (isodd(j) ? -staggered_field : staggered_field), "Sz", j
        end
    end
    H = MPO(ComplexF64, terms, sites)
    initial = MPS(ComplexF64, sites, [isodd(j) ? "Up" : "Dn" for j in 1:n])
    return (; sites, H, initial)
end
