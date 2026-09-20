#!/usr/bin/env julia
# Research-only analysis: no production backend import or Hamiltonian eigensolve.
const WORKER_START = time_ns()
using LinearAlgebra, SparseArrays, Serialization, SHA, TOML, Dates
include("../test/reference_ed.jl")
include("static18_transition_reference.jl")

const ROOT = normpath(joinpath(@__DIR__, ".."))
const PRIOR = "docs/research/data/p4_static18_validation.toml"
const PRIOR_HASH = "0d0e2c5b9ab23944babea2cd605ff5d712b23a0cb5988d5a544188bfae7ac702"
const PAYLOAD = "outputs/p4-static18-20260919-run2/reference_Q2.jls"
const PAYLOAD_HASH = "a5715b312c92eb4d5979d8da0599c5f7015d8905592b833f92f10588e570cdf7"
const SOURCES = ("examples/analyze_static18_neutral.jl", "examples/static18_transition_reference.jl",
    "examples/run_static18_neutral.py", "examples/run_static18.py", "examples/run_static27.py",
    "test/reference_ed.jl")
const LIMITS = Dict("residual"=>1e-10, "orthogonality"=>1e-10,
    "translation_closure"=>1e-9, "translation_unitarity"=>1e-9,
    "sum_rule"=>1e-10, "rotation_invariance"=>1e-10, "translation_covariance"=>1e-9,
    "parseval"=>1e-10, "spectral_weight_bound"=>1e-10, "calibration"=>1e-12)
file_hash(path) = bytes2hex(sha256(read(joinpath(ROOT,path))))
source_hashes() = Dict(p=>file_hash(p) for p in SOURCES)
elapsed() = (time_ns()-WORKER_START)/1e9
complex_record(A) = Dict("real"=>[real.(collect(row)) for row in eachrow(A)],
                         "imag"=>[imag.(collect(row)) for row in eachrow(A)])
canonical_bond(i,j,wy) = i<j ? (i,j,wy) : (j,i,-wy)

function persist(path, record)
    temporary, io = mktemp(dirname(path))
    try
        TOML.print(io,record;sorted=true)
        close(io)
        mv(temporary,path;force=true)
    finally
        isopen(io) && close(io)
        isfile(temporary) && rm(temporary)
    end
end

function load_input()
    file_hash(PRIOR)==PRIOR_HASH && file_hash(PAYLOAD)==PAYLOAD_HASH || error("input hash mismatch")
    prior = TOML.parsefile(joinpath(ROOT,PRIOR))
    prior["status"]=="passed_finite_cluster_comparison" || error("prior study did not pass")
    run = only(filter(r->r["Q"]==2,prior["runs"]))
    run["reference"]["payload"]["sha256"]==PAYLOAD_HASH || error("wrong payload association")
    run["reference_valid"] && run["reference"]["converged"] || error("unresolved prior reference")
    prior["extra_source_sha256"]["test/reference_ed.jl"]==file_hash("test/reference_ed.jl") ||
        error("independent reference changed since input generation")
    string(VERSION)==prior["runtime"]["julia"] || error("Julia serialization compatibility needs review")
    cfg = run["configuration"]
    cfg["Lx"]==2 && cfg["Ly"]==3 && cfg["N"]==18 && cfg["Q"]==2 &&
        cfg["ordering"]=="x_then_y_then_A_B_C" && cfg["axis_boundary"]=="open" &&
        cfg["circumference_boundary"]=="periodic" && cfg["gauge"]=="seam" &&
        cfg["charge_convention"]=="q=2Sz" && all(iszero,cfg["hz"]) || error("wrong configuration")
    sites,bonds = _reference_sites(2,3),reference_bonds(2,3)
    length(bonds)==30 && length(cfg["sites"])==18 || error("wrong geometry size")
    for (i,s) in enumerate(sites)
        old = cfg["sites"][i]
        (old["index"],old["x"],old["y"],old["sublattice"])==(i,s.x,s.y,string(s.sublattice)) &&
            norm(old["position"]-collect(s.position))<1e-14 || error("site geometry mismatch")
    end
    old_bonds = sort([canonical_bond(b["i"],b["j"],b["wy"]) for b in cfg["bonds"]])
    old_bonds==sort([(b.i,b.j,b.wy) for b in bonds]) &&
        all(b->b["family"]=="J1" && b["Jxy"]==b["Jz"]==1.0,cfg["bonds"]) ||
        error("NN bond geometry/couplings mismatch")
    # Deserialization is restricted to this hash-pinned, trusted local artifact.
    payload = deserialize(joinpath(ROOT,PAYLOAD))
    payload.Q==2 && payload.theta==0.0 && payload.values==run["reference"]["values"] ||
        error("payload metadata mismatch")
    payload.basis==_reference_basis(18;Q=2) && size(payload.vectors)==(43758,4) &&
        all(isfinite,payload.vectors) || error("invalid basis or vectors")
    return (;prior,run,payload,sites,bonds)
end

function translation_maps(basis,sites,bonds)
    site_map = [only(findall(t->t.x==s.x && t.y==mod(s.y+1,3) &&
        t.sublattice==s.sublattice,sites)) for s in sites]
    lookup = Dict(state=>i for (i,state) in enumerate(basis))
    state_map = [lookup[sum(((s>>(i-1))&UInt64(1))<<(site_map[i]-1) for i in 1:18)] for s in basis]
    pairs = Dict((b.i,b.j)=>i for (i,b) in enumerate(bonds))
    bond_map = [pairs[minmax(site_map[b.i],site_map[b.j])] for b in bonds]
    all(state_map[state_map[state_map[i]]]==i for i in eachindex(basis)) || error("T cubed is not identity")
    return (;site_map,state_map,bond_map)
end

function translation_subspace(V,p)
    TV = similar(V)
    TV[p,:] = V # Active T moves occupied sites y -> y+1; T|psi>=lambda|psi>.
    reduced = V'*TV
    return Dict("matrix"=>complex_record(reduced),"closure_residual"=>norm(TV-V*reduced),
        "unitarity_error"=>norm(reduced'*reduced-I),
        "eigenphase_over_pi"=>sort(angle.(eigvals(reduced))./pi),
        "cube_error"=>norm(reduced^3-I))
end

function orbit_fourier(A,p)
    seen = falses(length(p))
    records = Dict{String,Any}[]
    for i in eachindex(p)
        seen[i] && continue
        orbit = [i,p[i],p[p[i]]]
        length(unique(orbit))==3 && p[last(orbit)]==i || error("expected length-three operator orbit")
        seen[orbit] .= true
        # O_k=1/sqrt(3) sum_y exp(-ik*y) O_{T^y i}; no conjugation of coefficients.
        weights = [sum(abs2,[sum(cis(-2pi*m*y/3)*A[orbit[y+1],a] for y in 0:2)/sqrt(3)
            for a in axes(A,2)]) for m in (-1,0,1)]
        push!(records,Dict("operator_indices"=>orbit,"m"=>[-1,0,1],
            "k_over_pi"=>[-2/3,0.0,2/3],"projected_weight"=>weights))
    end
    totals = reduce(+,getindex.(records,"projected_weight"))
    return Dict("orbits"=>records,"m"=>[-1,0,1],"total_weight"=>totals,
        "parseval_error"=>abs(sum(totals)-sum(abs2,A)))
end

function transition_record(A,rotated,phased,p,ground_expectation,kind)
    C = conj(A)*transpose(A) # <g|O_i P_doublet O_j|g>; not full connected correlation.
    Cr = conj(rotated)*transpose(rotated)
    Cp = conj(phased)*transpose(phased)
    weights = real.(diag(C))
    means = real.(vec(ground_expectation))
    variance = kind=="sz" ? 0.25 .-means.^2 : 3/16 .-means./2 .-means.^2
    return Dict("amplitudes"=>complex_record(A),"projected_covariance"=>complex_record(C),
        "local_weight"=>weights,"total_weight"=>sum(weights),
        "ground_expectation"=>means,"ground_expectation_imaginary_error"=>maximum(abs,imag.(ground_expectation)),
        "local_connected_variance"=>variance,"local_variance_sum"=>sum(variance),
        "weight_minus_local_variance_max"=>maximum(weights-variance),
        "fraction_of_summed_local_variances"=>sum(weights)/sum(variance),
        "squared_singular_values"=>svdvals(A).^2,"sum_rule_error"=>maximum(abs,vec(sum(A;dims=1))),
        "rotation_covariance_error"=>norm(C-Cr),"ground_phase_covariance_error"=>norm(C-Cp),
        "translation_covariance_error"=>norm(C[p,p]-C),
        "fourier"=>orbit_fourier(A,p))
end

function main(output)
    isdir(output) && readdir(output)==["execution.toml"] || error("use launcher and fresh output")
    execution = TOML.parsefile(joinpath(output,"execution.toml"))
    execution["status"]=="running" && 0<execution["wall_limit_seconds"]<=120 || error("invalid supervisor")
    Threads.nthreads()==BLAS.get_num_threads()==1 || error("one Julia/BLAS thread required")
    before = source_hashes()
    record = Dict{String,Any}("schema_version"=>1,"status"=>"running",
        "recorded_at_utc"=>string(now(UTC)),"scope"=>"NN_N18_Q2_saved_neutral_doublet_transitions",
        "analysis_source_sha256"=>before,"analysis_git_revision"=>strip(read(`git -C $ROOT rev-parse HEAD`,String)),
        "analysis_worktree_dirty"=>!isempty(read(`git -C $ROOT status --porcelain`,String)),
        "runtime"=>Dict("julia"=>string(VERSION),"julia_threads"=>Threads.nthreads(),
            "blas_threads"=>BLAS.get_num_threads(),"LinearAlgebra"=>string(pkgversion(LinearAlgebra)),
            "Serialization"=>string(pkgversion(Serialization)),"SparseArrays"=>string(pkgversion(SparseArrays))),
        "limits"=>LIMITS,"wall_limit_seconds"=>execution["wall_limit_seconds"],
        "new_Hamiltonian_eigensolve"=>false,"new_DMRG"=>false,
        "small_matrix_decompositions"=>"1x1/2x2 projected translation eigenvalues and Nx2 transition SVD only",
        "bulk_gap"=>"not_established",
        "phase_identification"=>"not_attempted","execution_record"=>"execution.toml")
    function phase!(label)
        record["active_phase"] = label
        record["worker_elapsed_seconds"] = elapsed()
        persist(joinpath(output,"validation.toml"),record)
        println(label," elapsed=",round(elapsed();digits=2)); flush(stdout)
    end
    try
        phase!("load_and_check_provenance")
        input = load_input()
        ed,sites,bonds = input.payload,input.sites,input.bonds
        record["input"] = Dict("record_path"=>PRIOR,"record_sha256"=>PRIOR_HASH,
            "payload_path"=>PAYLOAD,"payload_sha256"=>PAYLOAD_HASH,
            "original_code"=>input.prior["code"],"original_runtime"=>input.prior["runtime"],
            "original_ed_solver"=>input.run["reference"],"configuration"=>input.run["configuration"],
            "interpretation"=>"original_ED_vectors_reanalyzed_not_relabelled_as_current_backend_states")
        phase!("independent_residual_translation_and_small_calibration")
        calibration = transition_reference_calibration()
        record["calibration"] = calibration
        calibration["passed"] || error("transition calibration failed")
        reference = reference_hamiltonian(2,3,0.0;Q=2,sparse_matrix=true)
        reference.basis==ed.basis || error("reference basis mismatch")
        V = ed.vectors
        residuals = [norm(reference.H*V[:,i]-ed.values[i]*V[:,i]) for i in 1:4]
        orthogonality = norm(V'*V-I)
        p = translation_maps(ed.basis,sites,bonds)
        h_error = norm(reference.H[p.state_map,p.state_map]-reference.H)
        ground_t = translation_subspace(V[:,1:1],p.state_map)
        doublet_t = translation_subspace(V[:,2:3],p.state_map)
        record["spectrum"] = Dict("values"=>ed.values,"residual_norms"=>residuals,
            "orthogonality_error"=>orthogonality,"neutral_spacing"=>ed.values[2]-ed.values[1],
            "doublet_splitting"=>ed.values[3]-ed.values[2],
            "next_returned_level_separation"=>ed.values[4]-ed.values[3],
            "completeness"=>"four_returned_levels_not_a_proof_of_full_low_spectrum",
            "basis_dimension"=>length(ed.basis),"Q"=>2,"theta"=>0.0)
        record["translation"] = Dict("active_convention"=>"site(x,y,s) maps to (x,(y+1) mod 3,s)",
            "state_convention"=>"T|psi> = exp(ik)|psi>","site_map"=>p.site_map,"bond_map"=>p.bond_map,
            "T_cubed_identity"=>true,"Hamiltonian_symmetry_error"=>h_error,
            "ground"=>ground_t,"doublet"=>doublet_t)
        observed_doublet = ed.values[2]-ed.values[1]>1e-8 && abs(ed.values[3]-ed.values[2])<1e-8 &&
            ed.values[4]-ed.values[3]>1e-8
        record["spectrum"]["cluster_tolerance"] = 1e-8
        record["spectrum"]["observed_doublet_separated"] = observed_doublet
        maximum(residuals)<=LIMITS["residual"] && orthogonality<=LIMITS["orthogonality"] &&
            observed_doublet && h_error==0 && all(t["closure_residual"]<=LIMITS["translation_closure"] &&
                max(t["unitarity_error"],t["cube_error"])<=LIMITS["translation_unitarity"] for t in (ground_t,doublet_t)) ||
            error("saved spectral/translation validation failed")
        phase!("transition_measurement")
        g,E = V[:,1],V[:,2:3]
        a = transition_amplitudes(g,E,ed.basis,18,bonds)
        means = transition_amplitudes(g,reshape(g,:,1),ed.basis,18,bonds)
        U = ComplexF64[1 im; im 1]/sqrt(2)
        rotated = transition_amplitudes(g,E*U,ed.basis,18,bonds)
        phased = transition_amplitudes(cis(0.37)*g,E,ed.basis,18,bonds)
        record["sz"] = transition_record(a.sz,rotated.sz,phased.sz,p.site_map,means.sz,"sz")
        record["bond"] = transition_record(a.bond,rotated.bond,phased.bond,p.bond_map,means.bond,"bond")
        spin = total_spin_diagnostics(V[:,1:3],ed.basis,18)
        record["total_spin_squared"] = spin.spin_squared
        record["spin_raising_norm_squared"] = spin.raising_norm_squared
        record["operators"] = Dict("sz"=>"physical Sz, hbar=1",
            "bond"=>"Si dot Sj, all J=1, each NN bond once; no gauge phase at theta=0",
            "amplitude"=>"A[i,alpha]=<e_alpha|O_i|g>, raw normalized input vectors",
            "covariance"=>"C[i,j]=sum_alpha conj(A[i,alpha])*A[j,alpha]=<g|O_i P_doublet O_j|g>",
            "fourier"=>"O_k=sum_y exp(-ik*y) O_(T^y i)/sqrt(3); relative momentum k_e-k_g=k",
            "sites"=>[Dict("index"=>i,"x"=>s.x,"y"=>s.y,"sublattice"=>string(s.sublattice),
                "position"=>collect(s.position)) for (i,s) in enumerate(sites)],
            "bonds"=>[Dict("index"=>i,"i"=>b.i,"j"=>b.j,"wy"=>b.wy) for (i,b) in enumerate(bonds)])
        checks = Dict("input_hash_and_geometry"=>true,"charge_and_complete_basis"=>true,
            "independent_residual"=>maximum(residuals)<=LIMITS["residual"],
            "orthogonality"=>orthogonality<=LIMITS["orthogonality"],
            "translation"=>true,"observed_doublet_separated"=>observed_doublet,
            "small_complex_operator_calibration"=>calibration["passed"])
        for kind in ("sz","bond")
            r = record[kind]
            checks[kind*"_sum_rule"] = r["sum_rule_error"]<=LIMITS["sum_rule"]
            checks[kind*"_rotation_invariance"] = max(r["rotation_covariance_error"],
                r["ground_phase_covariance_error"])<=LIMITS["rotation_invariance"]
            checks[kind*"_translation_covariance"] = r["translation_covariance_error"]<=LIMITS["translation_covariance"]
            checks[kind*"_parseval"] = r["fourier"]["parseval_error"]<=LIMITS["parseval"]
            checks[kind*"_spectral_weight_bound"] = r["weight_minus_local_variance_max"]<=LIMITS["spectral_weight_bound"] &&
                r["ground_expectation_imaginary_error"]<=LIMITS["spectral_weight_bound"] &&
                minimum(r["local_connected_variance"])>=-LIMITS["spectral_weight_bound"]
        end
        checks["analysis_sources_unchanged"] = before==source_hashes()
        checks["inputs_unchanged"] = file_hash(PRIOR)==PRIOR_HASH && file_hash(PAYLOAD)==PAYLOAD_HASH
        record["checks"] = checks
        record["status"] = all(values(checks)) ? "passed_saved_subspace_analysis" : "failed_analysis_checks"
    catch exception
        record["status"] = "failed"
        record["error_type"] = string(typeof(exception))
        rethrow()
    finally
        record["worker_elapsed_seconds"] = elapsed()
        persist(joinpath(output,"validation.toml"),record)
    end
    return record
end

if abspath(PROGRAM_FILE)==@__FILE__
    length(ARGS)==1 || error("expected fresh output directory")
    record = main(abspath(ARGS[1]))
    println("static18 neutral: ",record["status"])
    record["status"]=="passed_saved_subspace_analysis" || exit(1)
end
