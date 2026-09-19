#!/usr/bin/env julia
# Bounded interacting-cylinder accuracy study; no phase or plateau claim.
# julia --project=research --startup-file=no --threads=1 examples/validate_interacting18.jl
using SHA, TOML, Dates
const ROOT = normpath(joinpath(@__DIR__, ".."))
function active_environment_snapshot()
    project = Base.active_project()
    project === nothing && error("an explicit active project is required")
    project = abspath(project)
    manifest = Base.project_file_manifest_path(project)
    manifest !== nothing && isfile(project) && isfile(manifest) ||
        error("the active project needs an instantiated dependency manifest")
    manifest = abspath(manifest)
    return Dict{String,Any}("manifest"=>basename(manifest),
        "environment_sha256"=>Dict(
            "project:" * basename(project)=>bytes2hex(sha256(read(project))),
            "manifest:" * basename(manifest)=>bytes2hex(sha256(read(manifest)))))
end

const SOURCES = ["Project.toml", "research/Project.toml", "research/Manifest.toml",
    "research/Manifest-v1.13.toml",
    "src/KagomeDMRG.jl", "src/lattice.jl", "src/model.jl", "src/dmrg.jl",
    "src/observables.jl", "src/schmidt.jl", "src/checkpoint.jl", "src/continuation.jl",
    "test/reference_ed.jl", "test/reference_eigensolve.jl", "test/itensor_helpers.jl",
    "examples/validate_interacting18.jl"]
function source_hashes()
    hashes = Dict(path=>bytes2hex(sha256(read(joinpath(ROOT, path)))) for path in SOURCES)
    return merge(Dict("source_sha256"=>hashes), active_environment_snapshot())
end
const SOURCE_BEFORE = source_hashes()

using KagomeDMRG, ITensors, ITensorMPS, LinearAlgebra, Serialization
BLAS.set_num_threads(1)
ITensors.disable_threaded_blocksparse()
include(joinpath(ROOT,"test","reference_ed.jl"))
include(joinpath(ROOT,"test","reference_eigensolve.jl"))
include(joinpath(ROOT,"test","itensor_helpers.jl"))

const SOLVER = (; nsweeps=6, cutoff=0.0, noise=0.0, eigsolve_tol=1e-11,
    eigsolve_krylovdim=40, eigsolve_maxiter=20, measure_variance=true)
# These finite diagnostic limits are fixed before observing this study's paths.
# They do not select an expected pump value or lowest ED eigenstate.
const POLICY = FluxPolicy(min_overlap=0.90, max_density_change=0.05,
    max_entropy_change=0.10, max_schmidt_change=0.05, max_variance=1e-8,
    max_truncation_error=1e-10, max_sweep_energy_change=1e-9,
    max_cut_spread=1e-8, consistency_tol=1e-9)
const ED_LIMITS = Dict("energy_error_per_site"=>1e-8, "residual_norm"=>2e-6,
    "max_sz_error"=>1e-6, "ground_space_leakage"=>2e-5,
    "schmidt_density_error"=>1e-9, "variance_residual_squared_error"=>1e-9)

function run_study(directory; preflight=nothing)
    lattice = kagome_cylinder(2,3)
    record = Dict{String,Any}("schema_version"=>1,"status"=>"running",
        "recorded_at_utc"=>string(now(UTC)), "scope"=>"18_site_interacting_accuracy_study",
        "magnetization_plateau"=>"not_established", "phase_identification"=>"not_attempted",
        "quantization"=>"not_tested", "multiple_geometric_cuts"=>false,
        "bulk_region"=>"none_on_two_column_cylinder_all_sites_monitored",
        "configuration"=>KagomeDMRG._checkpoint_configuration(lattice,:seam,nothing),
        "runtime"=>merge(KagomeDMRG._checkpoint_runtime(),Dict(
            "julia_threads"=>Threads.nthreads(),"blas_threads"=>BLAS.get_num_threads())),
        "code"=>Dict{String,Any}("before"=>SOURCE_BEFORE,
            "provenance"=>KagomeDMRG._checkpoint_provenance()),
        "policy"=>KagomeDMRG._flux_policy_record(POLICY),
        "ed_limits"=>ED_LIMITS, "runs"=>Dict{String,Any}[],
        "references"=>Dict{String,Any}[], "comparisons"=>Dict{String,Any}())
    output = joinpath(directory,"validation.toml")
    persist() = KagomeDMRG._write_flux_record(output,record)
    persist()
    references = Dict{Float64,Any}()
    function reference(theta)
        return get!(references,theta) do
            # A different seed and block size from the independent audit.
            elapsed = @elapsed ed = reference_eigensystem(2,3,theta;seed=5678,blocksize=4)
            cluster = findall(x -> x-first(ed.values) < 1e-8,ed.values)
            resolved = length(cluster) < length(ed.values)
            info = Dict{String,Any}("theta"=>theta,"elapsed_seconds"=>elapsed,
                "values"=>ed.values,"residual_norms"=>ed.residual_norms,
                "converged"=>ed.converged,"orthogonality_error"=>ed.orthogonality_error,
                "observed_ground_cluster_size"=>length(cluster),
                "cluster_threshold"=>1e-8,"level_above_cluster_available"=>resolved,
                "settings"=>Dict(string(k)=>v for (k,v) in pairs(ed.settings)),
                "solver_info"=>Dict(string(k)=>v for (k,v) in pairs(ed.solver_info)))
            resolved && (info["finite_cluster_level_spacing"] = ed.values[length(cluster)+1]-first(ed.values))
            push!(record["references"],info)
            persist()
            ed.converged && resolved || error("independent ED reference is unresolved")
            # Trusted same-environment cache for independent validation; not
            # part of the portable TOML record or a production checkpoint.
            open(joinpath(directory,"reference_$(length(references)+1).jls"),"w") do io
                serialize(io,(;theta,basis=ed.basis,values=ed.values,vectors=ed.vectors,
                    residual_norms=ed.residual_norms,settings=ed.settings))
            end
            return ed
        end
    end
    function metrics(saved)
        ed = reference(saved.theta)
        state = sector_amplitudes(saved.psi,saved.sites,ed.basis)
        state ./= norm(state)
        cluster = findall(x -> x-first(ed.values) < 1e-8,ed.values)
        ground = ed.vectors[:,cluster]
        projected = ground*(ground'*state)
        projection_norm = norm(projected)
        matched_sz = reference_sz(projected,ed.basis,18)
        E = real(dot(state,ed.H*state))
        residual = norm(ed.H*state-E*state)
        diag = schmidt_diagnostics(saved.psi,9)
        transfer = only(spin_transfer(lattice,saved.baseline.sz,saved.diagnostics.sz))
        base_schmidt = schmidt_diagnostics(saved.baseline.psi,9)
        d = Dict{String,Any}("theta"=>saved.theta,"energy"=>saved.diagnostics.energy,
            "ed_energy"=>first(ed.values),
            "energy_error_per_site"=>abs(saved.diagnostics.energy-first(ed.values))/18,
            "mpo_reference_energy_difference"=>abs(saved.diagnostics.energy-E),
            "residual_norm"=>residual,"raw_variance"=>saved.diagnostics.variance,
            "variance_residual_squared_error"=>abs(saved.diagnostics.variance-residual^2),
            "ground_space_leakage"=>norm(state-projected),
            "ground_space_weight"=>projection_norm^2,
            "max_sz_error"=>maximum(abs.(saved.diagnostics.sz-matched_sz)),
            "maxlinkdim"=>maxlinkdim(saved.psi),"sz_profile"=>saved.diagnostics.sz,
            "schmidt_entropy"=>diag.entropy,"schmidt_mean_left_sz"=>diag.mean_left_sz,
            "schmidt_density_error"=>abs(diag.mean_left_sz-sum(saved.diagnostics.sz[1:9])),
            "transfer_left"=>transfer.left,"transfer_right"=>transfer.right,
            "transfer_total"=>transfer.total,
            "schmidt_transfer_right"=>-(diag.mean_left_sz-base_schmidt.mean_left_sz),
            "last_sweep_truncation_error"=>last(saved.diagnostics.max_truncation_errors))
        d["ed_comparison"] = all(d[key] <= limit for (key,limit) in ED_LIMITS) ? "passed" : "outside_limits"
        return d
    end
    starts = Dict{Tuple{Int,Int},Any}()
    function baseline(chi,seed)
        return get!(starts,(chi,seed)) do
            if preflight !== nothing && seed == 11
                root = dirname(abspath(preflight))
                old = TOML.parsefile(preflight)
                candidates = filter(r -> r["chi"]==chi && haskey(r,"checkpoint"),old["runs"])
                if length(candidates)==1
                    saved = load_checkpoint(joinpath(root,only(candidates)["checkpoint"]),lattice;
                        expected_theta=0.0)
                    expected = merge(SOLVER,(;maxdim=[chi],seed))
                    all(getproperty(saved.settings,k)==v for (k,v) in pairs(expected)) ||
                        error("preflight solver mismatch")
                    return KagomeDMRG._flux_result(saved)
                end
            end
            return run_dmrg(lattice,0.0;seed,maxdim=chi,SOLVER...)
        end
    end
    specs = [
        (;label="chi128_seed11_coarse",chi=128,seed=11,step=0.185,targets=[0.37,0.0]),
        (;label="chi512_seed11_coarse",chi=512,seed=11,step=0.185,targets=[0.37,0.0]),
        (;label="chi512_seed11_fine",chi=512,seed=11,step=0.0925,targets=[0.37,0.0]),
        (;label="chi512_seed29_coarse",chi=512,seed=29,step=0.185,targets=[0.37,0.0]),
        (;label="chi512_seed11_negative",chi=512,seed=11,step=0.185,targets=[-0.37,0.0])]
    try
        for spec in specs
            println("Starting ",spec.label)
            row = Dict{String,Any}("label"=>spec.label,"chi"=>spec.chi,"seed"=>spec.seed,
                "initial_step"=>spec.step,"targets"=>spec.targets,"status"=>"running",
                "points"=>Dict{String,Any}[])
            push!(record["runs"],row)
            persist()
            elapsed = @elapsed begin
                start = baseline(spec.chi,spec.seed)
                trajectory = continue_flux(lattice,spec.targets;start,output_root=directory,
                    policy=POLICY,initial_step=spec.step,min_step=spec.step/8,max_trials=24,
                    cuts=[1],diagnostic_bonds=[9],bulk_sites=1:18)
                row["status"] = string(trajectory.status)
                row["reason"] = trajectory.reason
                row["theta_path"] = trajectory.theta_path
                row["journal"] = relpath(trajectory.output_path,directory)
                journal = TOML.parsefile(trajectory.output_path)
                row["initial_diagnostics"] = journal["initial"]["diagnostics"]
                row["initial_reasons"] = journal["initial"]["reasons"]
                row["solver"] = journal["solver"]
                initial_path = joinpath(dirname(trajectory.output_path),journal["initial"]["checkpoint"])
                saved = load_checkpoint(initial_path,lattice;status=:trial)
                initial_metrics = metrics(saved)
                initial_metrics["point_status"] = journal["initial"]["status"]
                initial_metrics["checkpoint"] = relpath(initial_path,directory)
                push!(row["points"],initial_metrics)
                row["trials"] = [Dict("theta"=>t["theta"],"from_theta"=>t["from_theta"],
                    "status"=>t["status"],"reasons"=>t["reasons"],"diagnostics"=>t["diagnostics"])
                    for t in trajectory.trials]
                for trial in trajectory.trials
                    isempty(trial["checkpoint"]) && continue
                    path = joinpath(dirname(trajectory.output_path),trial["checkpoint"])
                    saved = load_checkpoint(path,lattice;status=:trial)
                    d = metrics(saved)
                    d["point_status"] = trial["status"]
                    d["checkpoint"] = relpath(path,directory)
                    push!(row["points"],d)
                end
            end
            row["elapsed_seconds"] = elapsed
            persist()
            println(spec.label," status=",row["status"]," points=",length(row["points"]),
                " elapsed_seconds=",elapsed)
        end
        byname = Dict(r["label"]=>r for r in record["runs"])
        # Compare common accepted fluxes without erasing direction/history.
        function compare(a,b)
            pairs = Dict{String,Any}[]
            for (index_a,p) in enumerate(a["points"])
                p["point_status"]=="accepted" || continue
                matches = [(j,q) for (j,q) in enumerate(b["points"]) if q["point_status"]=="accepted" &&
                    isapprox(q["theta"],p["theta"];atol=1e-13,rtol=0)]
                occurrence = count(q -> q["point_status"]=="accepted" &&
                    isapprox(q["theta"],p["theta"];atol=1e-13,rtol=0),a["points"][1:index_a])
                length(matches)>=occurrence || continue
                index_b,q = matches[occurrence]
                push!(pairs,Dict("theta"=>p["theta"],"point_index_a"=>index_a,
                    "point_index_b"=>index_b,"visit_number"=>occurrence,
                    "energy_difference"=>abs(p["energy"]-q["energy"]),
                    "max_sz_difference"=>maximum(abs.(p["sz_profile"]-q["sz_profile"])),
                    "right_transfer_difference"=>abs(p["transfer_right"]-q["transfer_right"]),
                    "schmidt_entropy_difference"=>abs(p["schmidt_entropy"]-q["schmidt_entropy"])))
            end
            return pairs
        end
        record["comparisons"]["step_halving"] = compare(byname["chi512_seed11_coarse"],byname["chi512_seed11_fine"])
        record["comparisons"]["seed_change"] = compare(byname["chi512_seed11_coarse"],byname["chi512_seed29_coarse"])
        positive = filter(p -> p["point_status"]=="accepted" && p["theta"]==0.37,
            byname["chi512_seed11_coarse"]["points"])
        negative = filter(p -> p["point_status"]=="accepted" && p["theta"]==-0.37,
            byname["chi512_seed11_negative"]["points"])
        if !isempty(positive) && !isempty(negative)
            p,n = first(positive),first(negative)
            sp = load_checkpoint(joinpath(directory,p["checkpoint"]),lattice;status=:trial)
            sn = load_checkpoint(joinpath(directory,n["checkpoint"]),lattice;status=:trial)
            vp = sector_amplitudes(sp.psi,sp.sites,reference(0.37).basis)
            vn = sector_amplitudes(sn.psi,sn.sites,reference(-0.37).basis)
            record["comparisons"]["opposite_flux_conjugacy"] = Dict(
                "energy_difference"=>abs(p["energy"]-n["energy"]),
                "max_sz_difference"=>maximum(abs.(p["sz_profile"]-n["sz_profile"])),
                "conjugate_overlap_amplitude"=>abs(dot(vn,conj(vp)))/(norm(vn)*norm(vp)))
        end
        for row in record["runs"]
            points = filter(p -> p["point_status"]=="accepted",row["points"])
            if length(points)>1 && last(points)["theta"]==0.0
                row["return_max_sz_difference"] = maximum(abs.(last(points)["sz_profile"]-first(points)["sz_profile"]))
                row["return_energy_difference"] = abs(last(points)["energy"]-first(points)["energy"])
                row["return_transfer_right"] = last(points)["transfer_right"]
            end
        end
        record["status"] = "completed_with_recorded_outcomes"
    catch exception
        exception isa InterruptException && rethrow()
        record["status"] = "error"
        record["exception_type"] = string(nameof(typeof(exception)))
        showerror(stderr,exception)
        println(stderr)
    end
    after = source_hashes()
    record["code"]["after"] = after
    record["code"]["status"] = after==SOURCE_BEFORE ? "verified" : "source_changed"
    after==SOURCE_BEFORE || (record["status"]="source_changed")
    persist()
    println("Saved interacting study: ",output)
    record["status"]=="completed_with_recorded_outcomes" || error("study incomplete; available outcomes retained")
end

length(ARGS)<=2 || error("usage: validate_interacting18.jl [new-output-directory] [trusted-preflight-record]")
directory = if isempty(ARGS)
    parent=joinpath(ROOT,"outputs")
    mkpath(parent)
    mktempdir(parent;prefix="p2-interacting18-",cleanup=false)
else
    path=abspath(ARGS[1])
    ispath(path) && error("output already exists; choose a new directory")
    mkpath(path)
    path
end
run_study(directory;preflight=length(ARGS)==2 ? ARGS[2] : nothing)
