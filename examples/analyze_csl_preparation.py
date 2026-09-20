#!/usr/bin/env python3
"""Compare the N36 Q0 chi256 zero-flux trial with one two-sweep preparation.

Usage: analyze_csl_preparation.py OUTPUT.json --old OLD.toml --new NEW.toml
       --execution NEW_EXECUTION.toml --config NEW_CONFIG.toml
       [--timed-out-variance]
Only saved measurements and hashes are checked; no Julia or MPS contraction runs.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import math
import os
from pathlib import Path
import sys
import tempfile
import tomllib

ROOT = Path(__file__).resolve().parent.parent
LIMITS = {"variance": 1e-3, "last_sweep_truncation_error": 1e-5,
          "sweep_energy_change": 1e-4}
INTEGRITY = {"checkpoint_overlap_error": 1e-10, "scaled_bond_energy_sum_error": 1e-9,
             "schmidt_density_error": 1e-9, "schmidt_probability_sum_error": 1e-10,
             "state_norm_error": 1e-10, "total_sz_error": 1e-9}
EXTRA_SOURCES = {"examples/research_csl.jl", "examples/run_research.py",
                 "examples/run_static18.py", "examples/run_static27.py"}


def require(condition, message):
    if not condition:
        raise ValueError(message)


def rootpath(value):
    path = Path(value)
    return (path if path.is_absolute() else ROOT / path).resolve()


def label(path):
    return os.path.relpath(path, ROOT)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def read_toml(path, evidence):
    payload = path.read_bytes()
    evidence[label(path)] = hashlib.sha256(payload).hexdigest()
    return tomllib.loads(payload.decode("utf-8"))


def number(value):
    require(type(value) in (int, float) and math.isfinite(value), "missing/nonfinite numeric measurement")
    return float(value)


def vector(values, size):
    require(isinstance(values, list) and len(values) == size, "profile length mismatch")
    return [number(value) for value in values]


def close(left, right, tolerance=1e-10):
    require(abs(number(left) - number(right)) <= tolerance, "saved/recomputed measurement mismatch")


def difference(old, new):
    require(len(old) == len(new) and bool(old), "incompatible profiles")
    delta = [number(b) - number(a) for a, b in zip(old, new)]
    return {"delta": delta, "max_abs": max(map(abs, delta)),
            "rms": math.sqrt(math.fsum(x * x for x in delta) / len(delta)),
            "mean": math.fsum(delta) / len(delta)}


def verify_snapshot(row, record, evidence):
    snapshot = row["snapshot"]
    path = rootpath(snapshot["path"])
    hashes = {name: sha256(path / name) for name in ("metadata.toml", "state.jls", "checksums.toml")}
    for name, key in (("metadata.toml", "metadata_sha256"), ("state.jls", "state_sha256"),
                      ("checksums.toml", "checksums_file_sha256")):
        require(hashes[name] == snapshot[key], "snapshot digest mismatch: " + name)
        evidence[label(path / name)] = hashes[name]
    checksums = read_toml(path / "checksums.toml", evidence)
    require(checksums["schema_version"] == 1 and set(checksums["files"]) == {"metadata.toml", "state.jls"},
            "unsupported checkpoint checksums")
    for name, check in checksums["files"].items():
        require(check["sha256"] == hashes[name] and check["bytes"] == (path / name).stat().st_size,
                "checkpoint checksum table mismatch")
    meta = read_toml(path / "metadata.toml", evidence)
    require(meta["format"] == "KagomeDMRG.local_checkpoint" and meta["schema_version"] == 1 and
            meta["status"] == "trial" and meta["theta_path"] == [0.0], "checkpoint is not a zero-flux trial")
    require(meta["configuration"] == record["configuration"] and
            meta["provenance"]["source_sha256"] == record["code"]["source_sha256"] and
            meta["runtime"] == {k: v for k, v in record["runtime"].items()
                                if k not in ("julia_threads", "blas_threads")} and
            meta["settings"] == row["settings"], "checkpoint identity mismatch")
    state = meta["state"]
    require(state["Q"] == 0 and state["theta"] == 0.0 and state["sz"] == row["sz_profile"] and
            state["sweep_energies"] == row["sweep_energies"] and
            state["max_truncation_errors"] == row["measured_truncation_errors"] and
            state["variance_measured"] is row["settings"]["measure_variance"], "checkpoint diagnostics mismatch")
    close(state["energy"], row["energy"], 1e-12)
    if state["variance_measured"]:
        close(state["variance"], row["raw_variance"], 1e-12)
    else:
        require("variance" not in state, "postprocessed variance was relabeled as a checkpoint measurement")
    return {"path": label(path), "sha256": hashes, "status": meta["status"]}


def validate_record(record, timed_out_variance=False):
    require(record["schema_version"] == 1 and record["N"] == 36 and record["Q"] == 0 and
            record["couplings"] == {"J1": 1.0, "J2": 0.5, "J3": 0.5}, "wrong CSL model/sector")
    if timed_out_variance:
        require(all(key not in record for key in ("configuration_and_source_unchanged", "parent_unchanged")),
                "timeout mode requires absent worker finalization flags")
    else:
        require(record["configuration_and_source_unchanged"] is True and record["parent_unchanged"] is True,
                "missing or failed source/configuration/parent finalization")
    require(record["integrity_limits"] == INTEGRITY, "integrity limits changed")
    cfg = record["configuration"]
    require(all(cfg[k] == v for k, v in {"Lx": 3, "Ly": 4, "N": 36, "Q": 0,
            "gauge": "seam", "ordering": "x_then_y_then_A_B_C", "charge_convention": "q=2Sz",
            "axis_boundary": "open", "circumference_boundary": "periodic"}.items()) and
            vector(cfg["hz"], 36) == [0.0] * 36, "wrong geometry/gauge/field")
    require(record["added_chirality_term"] is False and record["uniform_zeeman_field"] == 0.0 and
            record["flux_observations"] == [] and record["runtime"]["julia_threads"] == 1 and
            record["runtime"]["blas_threads"] == 1, "wrong field/flux/thread configuration")
    sources = record["code"]["source_sha256"]
    require(isinstance(sources, dict) and len(sources) == 180 and all(
        isinstance(h, str) and len(h) == 64 and all(c in "0123456789abcdef" for c in h)
        for h in sources.values()), "missing/incompatible backend source identities")
    require(set(record["extra_source_sha256"]) == EXTRA_SOURCES, "unexpected analysis source allowlist")


def measurements(row, record, timed_out_variance=False):
    if timed_out_variance:
        require(row["status"] == "running" and row["active_phase"] == "variance" and
                all(key not in row for key in ("raw_variance", "variance_source", "variance_roundoff_scale",
                    "second_moment_imaginary_part", "second_moment_imaginary_part_measured",
                    "integrity_passed", "integrity_failures", "accuracy_status", "change_from_parent")),
                "not the supported incomplete variance-phase outcome")
        require(set(row["measurement_seconds"]) == {"schmidt", "bond_energies", "scalar_chirality"} and
                all(number(value) >= 0 for value in row["measurement_seconds"].values()),
                "pre-variance measurements did not complete")
    else:
        require(row["integrity_passed"] is True and row["integrity_failures"] == [] and
                row["active_phase"] == "diagnostics_completed", "incomplete or unverified diagnostics")
    require(row["checkpoint_overlap_measured"] is True, "checkpoint reload overlap was not measured")
    for key, limit in INTEGRITY.items():
        require(0 <= number(row[key]) <= limit, "integrity limit failed: " + key)
    require(row["theta"] == 0.0 and row["maxlinkdim"] == 256, "wrong theta/actual chi")
    settings = row["settings"]
    require(settings["nsweeps"] == 2 and settings["maxdim"] == [256] and
            settings["cutoff"] == settings["noise"] == 0 and
            settings["initialization"] == "provided_mps", "unexpected solver settings")
    sz = vector(row["sz_profile"], 36)
    bonds = vector(row["bond_energy_profile"], len(record["configuration"]["bonds"]))
    chirality = vector(row["chirality_profile"], len(record["triangles"]))
    energy = number(row["energy"])
    variance = None if timed_out_variance else number(row["raw_variance"])
    truncation = vector(row["measured_truncation_errors"], 2)
    sweeps = vector(row["sweep_energies"], 2)
    require(all(0 <= t <= 1 for t in truncation), "invalid truncation measurement")
    close(math.fsum(sz), 0, 1e-9)
    close(math.fsum(bonds), energy, 1e-9 * max(1, abs(energy)))
    close(math.fsum(chirality) / len(chirality), row["chirality_mean"])
    columns = [math.fsum(sz[i:i + 12]) for i in range(0, 36, 12)]
    for actual, saved in zip(columns, vector(row["column_sz"], 3)):
        close(actual, saved)
    sweep_change = max(abs(sweeps[1] - sweeps[0]), abs(energy - sweeps[-1]))
    close(sweep_change, row["sweep_energy_change"], 1e-12)
    roundoff = measured_imag = None
    if not timed_out_variance:
        roundoff = 100 * sys.float_info.epsilon * max(1, energy * energy)
        close(roundoff, row["variance_roundoff_scale"], 1e-20)
        require(variance >= -roundoff, "negative variance outside roundoff")
        measured_imag = row["second_moment_imaginary_part_measured"]
        require(type(measured_imag) is bool, "missing variance measurement provenance")
        if measured_imag:
            require(0 <= number(row["second_moment_imaginary_part"]) <= roundoff, "nonreal HdaggerH moment")
        else:
            require("second_moment_imaginary_part" not in row, "absent imaginary part relabeled as measured")
    schmidt = []
    require([s["bond"] for s in row["schmidt"]] == [12, 24], "wrong Schmidt cuts")
    for entry in row["schmidt"]:
        probabilities = vector(entry["probabilities"], len(entry["left_q"]))
        charges = vector(entry["left_q"], len(probabilities))
        require(bool(probabilities) and all(p >= 0 for p in probabilities) and
                all(q.is_integer() for q in charges), "invalid Schmidt probabilities/charges")
        close(math.fsum(probabilities), 1, 1e-10)
        mean_sz = math.fsum(p * q / 2 for p, q in zip(probabilities, charges))
        entropy = -math.fsum(p * math.log(p) for p in probabilities if p > 0)
        close(mean_sz, entry["mean_left_sz"], 1e-9)
        close(mean_sz, math.fsum(sz[:entry["bond"]]), 1e-9)
        close(entropy, entry["entropy"], 1e-9)
        schmidt.append({"bond": entry["bond"], "entropy": number(entry["entropy"]),
                        "mean_left_sz": number(entry["mean_left_sz"]), "mean_left_q": 2 * mean_sz,
                        "variance_left_sz": number(entry["variance_left_sz"])})
    return {"energy": energy, "variance": variance, "last_sweep_truncation_error": truncation[-1],
            "sweep_energy_change": sweep_change, "sweep_energies": sweeps,
            "measured_truncation_errors": truncation, "sz_profile": sz, "column_sz": columns,
            "bond_energy_profile": bonds, "chirality_profile": chirality,
            "chirality_mean": number(row["chirality_mean"]), "schmidt": schmidt,
            "variance_measurement": {"solver_measure_variance": settings["measure_variance"],
                "source": row.get("variance_source"), "roundoff_scale": roundoff,
                "status": "missing_timeout_during_postcheckpoint_HdaggerH" if timed_out_variance else "measured",
                "second_moment_imaginary_part_measured": measured_imag,
                "second_moment_imaginary_part": row.get("second_moment_imaginary_part")}}


def analyze(old_path, new_path, execution_path, config_path, *, timed_out_variance=False):
    evidence = {}
    require(old_path != new_path, "old and new validations must be distinct")
    old, new = (read_toml(p, evidence) for p in (old_path, new_path))
    validate_record(old)
    validate_record(new, timed_out_variance)
    require(old["mode"] == "flux" and old["status"] == "unresolved_flux" and
            old["continuation"]["status"] == "unresolved" and
            old["continuation"]["reason"] == "initial_diagnostics" and
            old["continuation"]["theta_path"] == [0.0] and
            old["continuation"]["accepted_checkpoints"] == [] and
            old["continuation"]["last_accepted_checkpoint"] == "", "wrong old flux-start outcome")
    require(new["mode"] == "prepare" and len(new["preparation_batches"]) == 1,
            "new preparation must contain exactly one batch")
    if timed_out_variance:
        require(new["status"] == "running" and "flux_readiness" not in new and
                "latest_trial_checkpoint" not in new, "worker finalized outside the supported timeout phase")
    else:
        require(new["status"] == "completed_preparation_accuracy_unestablished" and
                new["flux_readiness"] == "not_certified_by_preparation_pilot", "new preparation has not completed")
    before, after = old["flux_start"], new["preparation_batches"][0]
    require(before["status"] == "see_core_initial_policy_decision" and before["theta_path"] == [0.0] and
            after["status"] == ("running" if timed_out_variance else "completed_accuracy_unestablished") and after["batch"] == 1 and
            after["sweeps_in_this_invocation"] == 2, "wrong selected rows/sweep count")
    for key in ("configuration", "runtime", "couplings", "triangles", "cuts", "schmidt_bonds", "center_sites"):
        require(old[key] == new[key], "old/new mismatch: " + key)
    require(old["code"]["source_sha256"] == new["code"]["source_sha256"], "backend source mismatch")
    require(old["extra_source_sha256"]["examples/research_csl.jl"] ==
            new["extra_source_sha256"]["examples/research_csl.jl"], "measurement worker source mismatch")
    require(old["config"]["solver"] == new["config"]["solver"] and
            {k: v for k, v in before["settings"].items() if k != "measure_variance"} ==
            {k: v for k, v in after["settings"].items() if k != "measure_variance"}, "solver mismatch")
    require(before["settings"]["measure_variance"] is True and
            before["variance_source"] == "run_dmrg_measured_variance" and
            before["second_moment_imaginary_part_measured"] is False and
            after["settings"]["measure_variance"] is False, "wrong variance measurement settings")
    if not timed_out_variance:
        require(after["variance_source"] == "separate_HdaggerH_contraction_after_trial_checkpoint" and
                after["second_moment_imaginary_part_measured"] is True, "wrong new variance measurement path")
    cfg = read_toml(config_path, evidence)
    require(cfg == new["config"] and evidence[label(config_path)] == new["config_sha256"] and
            cfg["preparation"] == {"batches": 1, "batch_sweeps": 2} and
            cfg["start_status"] == "trial", "explicit config does not match new run")
    snapshots = {"old": verify_snapshot(before, old, evidence), "new": verify_snapshot(after, new, evidence)}
    parent = new["parent"]
    for key in ("path", "metadata_sha256", "state_sha256", "checksums_file_sha256"):
        require(parent[key] == before["snapshot"][key], "new parent does not match old selected checkpoint")
    require(parent["status"] == "trial" and parent["theta_path"] == [0.0] and
            parent["solver"] == before["settings"] and parent["variance_measured"] is True and
            rootpath(cfg["start_checkpoint"]) == rootpath(parent["path"]) == rootpath(after["parent_checkpoint"]),
            "new preparation parent linkage mismatch")
    close(parent["energy"], before["energy"], 1e-12)
    require(snapshots["old"]["path"] != snapshots["new"]["path"] and
            (new_path.parent / after["checkpoint"]).resolve() == rootpath(snapshots["new"]["path"]),
            "wrong new checkpoint selection")
    if not timed_out_variance:
        require(new["latest_trial_checkpoint"] == after["checkpoint"], "wrong latest completed trial")
    execution = read_toml(execution_path, evidence)
    require(execution_path == new_path.parent / new["execution_record"] and execution["schema_version"] == 1 and
            execution["worker_exit_confirmed"] is True and execution["julia_threads"] == execution["blas_threads"] == 1 and
            execution["termination_grace_included_in_wall_limit"] is True and
            0 < number(execution["wall_elapsed_seconds"]) <= number(execution["wall_limit_seconds"]) <= 600 and
            execution["wall_limit_seconds"] == new["wall_limit_seconds"], "new execution termination/bound mismatch")
    if timed_out_variance:
        require(execution["status"] == "timed_out" and execution["worker_exit_code"] in (-15, -9) and
                execution["wall_limit_seconds"] - execution["termination_grace_seconds"] <=
                number(execution["sigterm_elapsed_seconds"]) <= number(execution["sigkill_elapsed_seconds"]) <=
                execution["wall_elapsed_seconds"] and
                0 < number(new["worker_elapsed_seconds"]) < execution["sigterm_elapsed_seconds"],
                "not an externally terminated variance timeout")
    else:
        require(execution["status"] == "exited" and execution["worker_exit_code"] == 0,
                "new execution did not complete successfully")
    command = execution["command"]
    require(len(command) == 7 and command[:4] == ["julia", "--project=.", "--startup-file=no", "--threads=1"] and
            rootpath(command[4]) == ROOT / "examples/research_csl.jl" and
            rootpath(command[5]) == new_path.parent and rootpath(command[6]) == config_path,
            "execution command is not bound to the explicit new output/config")
    archive = new_path.parent / "analysis-sources"
    require(sha256(archive / "config.toml") == new["config_sha256"], "archived config hash mismatch")
    evidence[label(archive / "config.toml")] = new["config_sha256"]
    for relative, expected in new["extra_source_sha256"].items():
        require(sha256(archive / relative) == expected, "archived analysis source mismatch: " + relative)
        evidence[label(archive / relative)] = expected
    current_extra_matches = {relative: sha256(ROOT / relative) == expected
                             for relative, expected in new["extra_source_sha256"].items()}
    a, b = measurements(before, old), measurements(after, new, timed_out_variance)
    parent_change = {"energy": b["energy"] - a["energy"], "max_center_sz_change":
                     max(abs(b["sz_profile"][i] - a["sz_profile"][i]) for i in range(12, 24)),
                     "overlap": None if timed_out_variance else number(after["change_from_parent"]["overlap"])}
    if not timed_out_variance:
        close(after["change_from_parent"]["energy"], parent_change["energy"], 1e-12)
        close(after["change_from_parent"]["max_center_sz_change"], parent_change["max_center_sz_change"], 1e-12)
        require(0 <= parent_change["overlap"] <= 1 + 1e-10, "invalid parent overlap")
    policy = old["config"]["policy"]
    require({"variance": policy["max_variance"], "last_sweep_truncation_error": policy["max_truncation_error"],
             "sweep_energy_change": policy["max_sweep_energy_change"]} == LIMITS, "old readiness limits changed")
    policy_checks = {key: {"limit": limit, "old_value": a[key], "new_value": b[key],
                          "old_passed": abs(a[key]) <= limit,
                          "new_passed": None if b[key] is None else abs(b[key]) <= limit}
                     for key, limit in LIMITS.items()}
    decisions = [item["new_passed"] for item in policy_checks.values()]
    combined = False if False in decisions else None if None in decisions else True
    schmidt_changes = [{"bond": s["bond"], **{key + "_change": t[key] - s[key]
                        for key in ("entropy", "mean_left_sz", "mean_left_q", "variance_left_sz")}}
                       for s, t in zip(a["schmidt"], b["schmidt"])]
    evidence[label(Path(__file__).resolve())] = sha256(Path(__file__).resolve())
    for relative, expected in evidence.items():
        require(sha256(rootpath(relative)) == expected, "input changed during comparison: " + relative)
    return {"schema_version": 1, "recorded_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
            "analysis_mode": "timed_out_variance" if timed_out_variance else "completed_preparation",
            "study": "N36 Q0 J1=1 J2=J3=0.5 chi256 zero-flux preparation; two additional sweeps",
            "input_and_analysis_sha256": evidence, "old_validation": label(old_path),
            "new_validation": label(new_path), "new_execution": label(execution_path), "new_config": label(config_path),
            "verification_scope": "saved measurements, metadata and checkpoint payload hashes; no deserialization or new contraction; backend source digest maps compared, backend files not rehashed here",
            "configuration": new["configuration"], "triangles": new["triangles"],
            "backend_source_sha256": new["code"]["source_sha256"], "runtime": new["runtime"],
            "extra_source_sha256": {"old": old["extra_source_sha256"], "new": new["extra_source_sha256"]},
            "settings": {"old": before["settings"], "new": after["settings"]}, "checkpoints": snapshots,
            "source_and_parent_finalization_verified": None if timed_out_variance else True,
            "worker_recorded_flags": {"old": {key: old.get(key) for key in
                ("configuration_and_source_unchanged", "parent_unchanged")},
                "new": {key: new.get(key) for key in ("configuration_and_source_unchanged", "parent_unchanged")},
                "new_row": {key: after.get(key) for key in ("integrity_passed", "integrity_failures", "accuracy_status")}},
            "external_post_run_audit": {"parent_and_new_checkpoint_hashes_match": True,
                "archived_analysis_sources_match": True, "current_extra_source_hashes_match": current_extra_matches,
                "available_saved_integrity_checks_passed": True,
                "qualification": "Python checks after worker exit; these do not supply missing worker finalization flags or establish that source files were unchanged throughout execution"},
            "preserved_status": {"old_worker": old["status"], "old_row": before["status"],
                "old_continuation": old["continuation"], "new_worker": new["status"], "new_row": after["status"],
                "new_flux_readiness": new.get("flux_readiness"), "new_active_phase": after["active_phase"],
                "execution": execution["status"], "worker_exit_code": execution["worker_exit_code"]},
            "sweeps": {"additional_sweeps_verified_from_new_record": 2, "old_selected_batch_sweeps": 2,
                "absolute_cumulative_sweeps_verified_from_these_inputs": None,
                "historical_chain_context": {"old": 8, "new": 10,
                    "basis": "earlier prepare 4 + refine 2 + old flux-start 2, then new prepare 2; earlier chain not revalidated by this CLI"}},
            "old": a, "new": b, "change_convention": "new minus old; same site/bond/triangle ordering",
            "scalar_changes": {key: None if b[key] is None else b[key] - a[key] for key in
                               ("energy", "variance", "last_sweep_truncation_error", "sweep_energy_change", "chirality_mean")},
            "profile_changes": {key: difference(a[key], b[key]) for key in
                                ("sz_profile", "bond_energy_profile", "chirality_profile", "column_sz")},
            "schmidt_changes": schmidt_changes, "change_from_parent_recorded": after.get("change_from_parent"),
            "change_from_parent_external_comparison": parent_change,
            "checkpoint_reload_overlap_error": after["checkpoint_overlap_error"],
            "units": {"energy": "J", "variance": "J^2", "Sz": "hbar=1", "left_q": "2*Sz", "entropy": "natural logarithm"},
            "old_policy": policy, "readiness_criteria": policy_checks,
            "three_scalar_criteria_passed": combined,
            "readiness_status": "not_eligible_on_available_criteria" if combined is False else "not_certified",
            "actual_flux_policy_acceptance": False,
            "claims": {"accuracy": "unestablished", "matched_accuracy_inferred": False,
                "bulk_convergence": "not_established", "phase_identification": "not_attempted",
                "nonzero_flux_attempted": False, "nonzero_pump_calibrated": False,
                "readiness_comparison": "three scalar checks only; prepare does not execute FluxPolicy or promote a trial to accepted",
                "variance": "Hamiltonian dispersion, not a ground-state energy error bound"},
            "new_resources": {key: execution.get(key) for key in ("wall_limit_seconds", "wall_elapsed_seconds",
                "sigterm_elapsed_seconds", "sigkill_elapsed_seconds",
                "worker_cpu_total_seconds", "worker_peak_rss_bytes", "resource_measurement_status")}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    for name in ("old", "new", "execution", "config"):
        parser.add_argument("--" + name, type=Path, required=True)
    parser.add_argument("--timed-out-variance", action="store_true",
                        help="allow only a saved two-sweep trial terminated during its variance contraction")
    args = parser.parse_args()
    temporary = None
    try:
        output = args.output.resolve()
        require(output.suffix == ".json" and not output.exists(), "use a new .json output path")
        result = analyze(*(getattr(args, name).resolve() for name in ("old", "new", "execution", "config")),
                         timed_out_variance=args.timed_out_variance)
        encoded = json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
        output.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=output.parent, delete=False) as stream:
            temporary = Path(stream.name)
            stream.write(encoded)
            stream.flush()
            os.fsync(stream.fileno())
        os.link(temporary, output)  # Atomic publication; unlike replace, never overwrites another result.
        print(label(output))
        return 0
    except (OSError, KeyError, TypeError, ValueError) as error:
        print(f"CSL comparison: {error}", file=sys.stderr)
        return 2
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)


if __name__ == "__main__":
    sys.exit(main())
