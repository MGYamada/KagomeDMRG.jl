#!/usr/bin/env python3
"""Compare NN N27 Q5 chi256 cumulative8-to10, with Q1/Q3 cumulative8 fixed.

Only finalized, verified saved trials enter this comparison. No MPS is loaded
or contracted, and numerical variation is not a rigorous boundary error bar.
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

sys.dont_write_bytecode = True
import analyze_sector_refinement as refinement
import summarize_static_campaign as campaign

LIMITS = refinement.LIMITS
ROLES = {"fixed_q1": (1, 8), "fixed_q3": (3, 8), "old_q5": (5, 8), "new_q5": (5, 10)}
COMPLETE = "completed_bounded_case_accuracy_separate"
require, digest, label, rootpath = refinement.require, refinement.digest, refinement.label, refinement.rootpath


def load_input(path: Path, role: str) -> dict:
    path = path.resolve()
    payload = path.read_bytes()
    record = tomllib.loads(payload.decode("utf-8"))
    require(record.get("status") == COMPLETE, f"{role}: worker is not finalized and complete")
    q, sweeps = ROLES[role]
    require(record.get("schema_version") == 1 and record.get("N") == 27 and record.get("Q") == q,
            f"{role}: wrong schema/N/Q")
    require(record.get("sources_unchanged") is True and record.get("parent_unchanged") is True,
            f"{role}: source/parent finalization absent or failed")
    # This one historical record predates the backend_unchanged field. Preserve
    # its absence; matching saved backend maps does not recreate that check.
    backend_flag = record.get("backend_unchanged")
    require(backend_flag is True or (role == "fixed_q3" and "backend_unchanged" not in record),
            f"{role}: backend finalization absent or failed")
    geometry = campaign._geometry(record["configuration"], record)
    require((geometry["Lx"], geometry["Ly"]) == (3, 3) and geometry.get("gauge") == "seam"
            and geometry.get("ordering") == "x_then_y_then_A_B_C", f"{role}: wrong geometry/conventions")
    cfg, solver = record["config"], record["solver"]
    require(all(cfg[k] == v for k, v in {"Lx": 3, "Ly": 3, "Q": q, "maxdim": 256}.items())
            and cfg["stationarity"] == LIMITS, f"{role}: configuration or limits changed")
    require(solver["maxdim"] == 256 and solver["nsweeps"] == 2
            and solver["cutoff"] == solver["noise"] == 0, f"{role}: wrong solver")
    rows = record["batches"]
    require(rows and all(campaign._eligible_batch(row) for row in rows), f"{role}: incomplete batch")
    chosen = [row for row in rows if row["cumulative_sweeps"] == sweeps]
    require(len(chosen) == 1 and max(row["cumulative_sweeps"] for row in rows) == sweeps
            and record["completed_sweeps"] == sweeps, f"{role}: required final cumulative{sweeps} row absent")
    row = chosen[0]
    require(row["maxlinkdim"] == 256, f"{role}: actual chi is not 256")
    energy, variance = (campaign._number(row[k], k) for k in ("energy", "variance"))
    sz = campaign._vector(row["sz_profile"], 27, "sz_profile")
    bonds = campaign._vector(row["bond_energy"], len(geometry["bonds"]), "bond_energy")
    truncation = campaign._vector(row["measured_truncation_errors"], 2, "truncation")
    sweep_energies = campaign._vector(row["sweep_energies"], 2, "sweep_energies")
    require(all(0 <= x <= 1 for x in truncation) and abs(math.fsum(sz) - q / 2) <= 1e-10
            and abs(math.fsum(bonds) - energy) <= 1e-8, f"{role}: observable sum rule failed")
    floor = 100 * sys.float_info.epsilon * max(1, energy * energy)
    require(variance >= -floor, f"{role}: negative variance beyond roundoff")
    cp = rootpath(row["checkpoint"])
    hashes = {name: digest(cp / name) for name in ("metadata.toml", "state.jls", "checksums.toml")}
    require(hashes["metadata.toml"] == row["checkpoint_metadata_sha256"]
            and hashes["state.jls"] == row["checkpoint_payload_sha256"], f"{role}: checkpoint digest mismatch")
    checksums = tomllib.loads((cp / "checksums.toml").read_text())
    for name in ("metadata.toml", "state.jls"):
        require(checksums["files"][name] == {"sha256": hashes[name], "bytes": (cp / name).stat().st_size},
                f"{role}: checkpoint checksum table mismatch")
    meta = tomllib.loads((cp / "metadata.toml").read_text())
    require(meta["status"] == "trial" and meta["theta_path"] == [0.0]
            and meta["configuration"] == record["configuration"], f"{role}: checkpoint configuration mismatch")
    state = meta["state"]
    require(state["Q"] == q and state["theta"] == 0 and state["energy"] == energy and state["sz"] == sz
            and state["sweep_energies"] == sweep_energies and state["max_truncation_errors"] == truncation,
            f"{role}: checkpoint measurements mismatch")
    backend = record["code"]["source_sha256"]
    require(isinstance(backend, dict) and bool(backend) and meta["provenance"]["source_sha256"] == backend,
            f"{role}: checkpoint/backend identity mismatch")
    require(meta["runtime"] == {k: v for k, v in record["runtime"].items()
                                if k not in ("julia_threads", "blas_threads")}, f"{role}: checkpoint runtime mismatch")
    summary = {"validation_path": label(path), "validation_sha256": hashlib.sha256(payload).hexdigest(),
               "case_id": record["case_id"], "Q": q, "N": 27, "cumulative_sweeps": sweeps,
               "configured_maxdim": 256, "actual_maxlinkdim": row["maxlinkdim"],
               "selected_batch": row["batch"], "energy": energy, "variance": variance,
               "variance_per_site": variance / 27, "last_sweep_truncation_error": truncation[-1],
               "sz_profile": sz, "bond_energy": bonds, "sweep_energies": sweep_energies,
               "checkpoint": label(cp), "checkpoint_sha256": hashes,
               "source_finalization": {k: record.get(k) for k in
                                       ("sources_unchanged", "backend_unchanged", "parent_unchanged")},
               "quality": campaign._quality_summary(row), "backend_source_sha256": backend,
               "recorded_analysis_source_sha256": record["analysis_source_sha256"]}
    return {"path": path, "record": record, "row": row, "geometry": geometry, "summary": summary}


def check_new_run(old: dict, new: dict, execution_path: Path) -> dict:
    record, before = new["record"], old["summary"]
    cfg, parent = record["config"], record["parent"]
    require(cfg["initialization"] == "resume" and cfg["initial_completed_sweeps"] == 8
            and parent["initial_completed_sweeps"] == 8 and cfg["batches"] == 1 and cfg["batch_sweeps"] == 2,
            "new Q5 must be one two-sweep batch resumed from cumulative8")
    require(parent.get("strict_source_runtime_configuration_load") is True, "strict parent load absent")
    cp, evidence = rootpath(before["checkpoint"]), rootpath(parent["record"])
    require(rootpath(cfg["parent_checkpoint"]) == rootpath(parent["checkpoint"]) == cp,
            "parent checkpoint path mismatch")
    require(rootpath(cfg["parent_record"]) == evidence and digest(evidence) == before["validation_sha256"]
            == parent["record_sha256"] == cfg["parent_record_sha256"], "parent evidence digest mismatch")
    require(cfg["parent_metadata_sha256"] == before["checkpoint_sha256"]["metadata.toml"],
            "parent metadata digest mismatch")
    expected = {cp / name: sha for name, sha in before["checkpoint_sha256"].items()}
    expected[evidence] = before["validation_sha256"]
    require({rootpath(k): v for k, v in parent["checkpoint_sha256"].items()} == expected,
            "parent file digest table mismatch")
    require(abs(parent["loaded_energy"] - before["energy"]) <= 1e-12, "loaded parent energy mismatch")
    execution_path = execution_path.resolve()
    require(execution_path == campaign._companion_paths(new["path"])[0], "execution is not the new validation companion")
    execution = tomllib.loads(execution_path.read_text())
    require(execution.get("status") == "exited" and execution.get("worker_exit_confirmed") is True
            and execution.get("worker_exit_code") == 0, "new worker has not exited successfully")
    require(0 < execution["wall_limit_seconds"] <= 600 and execution["wall_elapsed_seconds"] <= execution["wall_limit_seconds"]
            and record["wall_limit_seconds"] == execution["wall_limit_seconds"], "new worker exceeded its wall budget")
    require(execution["julia_threads"] == execution["blas_threads"] == 1
            and record["runtime"]["julia_threads"] == record["runtime"]["blas_threads"] == 1, "wrong thread allocation")
    # The checkpoint identifies the original run directory even if validation
    # and execution were later copied byte-for-byte into docs/research/data.
    run_root = rootpath(new["summary"]["checkpoint"]).parent.parent
    require(digest(run_root / "validation.toml") == new["summary"]["validation_sha256"]
            and digest(run_root / "execution.toml") == digest(execution_path), "new run/archive binding mismatch")
    archived = run_root / "analysis-sources"
    require(digest(archived / "config.toml") == record["config_sha256"]
            and tomllib.loads((archived / "config.toml").read_text()) == cfg, "archived configuration mismatch")
    require(all(digest(archived / name) == sha for name, sha in record["analysis_source_sha256"].items()),
            "archived worker sources mismatch")
    return {"parent_link": "verified_cumulative8_to10", "parent_validation_sha256": before["validation_sha256"],
            "execution_path": label(execution_path), "execution_sha256": digest(execution_path),
            "execution": campaign._execution_summary(execution), "archived_config_and_analysis_sources_verified": True}


def analyze(paths: dict[str, Path], execution_path: Path) -> dict:
    # Read the new point first, so a live/incomplete worker is rejected before
    # any checkpoint hashing or analysis of prior trials.
    cases = {"new_q5": load_input(paths["new_q5"], "new_q5")}
    cases.update({role: load_input(paths[role], role) for role in ROLES if role != "new_q5"})
    reference = cases["fixed_q3"]
    geometry = {k: v for k, v in reference["geometry"].items() if k != "Q"}
    for role, case in cases.items():
        require({k: v for k, v in case["geometry"].items() if k != "Q"} == geometry, f"{role}: geometry mismatch")
        require(case["record"]["solver"] == reference["record"]["solver"], f"{role}: solver mismatch")
        require(case["record"]["runtime"] == reference["record"]["runtime"], f"{role}: runtime mismatch")
        require(case["summary"]["backend_source_sha256"] == reference["summary"]["backend_source_sha256"],
                f"{role}: backend source mismatch")
    verified = check_new_run(cases["old_q5"], cases["new_q5"], execution_path)
    rows = {role: case["summary"] for role, case in cases.items()}
    old, new, raw = rows["old_q5"], rows["new_q5"], cases["new_q5"]["row"]
    delta_e = new["energy"] - old["energy"]
    optimizer_error = abs(new["energy"] - new["sweep_energies"][-1])
    last_sweep = max(abs(new["sweep_energies"][-1] - new["sweep_energies"][-2]), optimizer_error)
    values = {"energy_per_site": max(abs(delta_e), last_sweep) / 27,
              "sz_profile": max(abs(x - y) for x, y in zip(old["sz_profile"], new["sz_profile"])),
              "bond_profile": max(abs(x - y) for x, y in zip(old["bond_energy"], new["bond_energy"])),
              "variance_per_site": abs(new["variance"]) / 27, "truncation": new["last_sweep_truncation_error"]}
    criteria = {key: {"value": value, "limit": LIMITS[key], "passed": value <= LIMITS[key]}
                for key, value in values.items()}
    passed = all(item["passed"] for item in criteria.values())
    for key, value in {"energy_change_from_previous": delta_e, "last_optimizer_energy_error": optimizer_error,
                       "last_sweep_energy_change": last_sweep, "max_sz_change": values["sz_profile"],
                       "max_bond_change": values["bond_profile"], "variance_per_site": new["variance"] / 27}.items():
        require(abs(campaign._number(raw[key], key) - value) <= 5e-13, f"worker arithmetic mismatch: {key}")
    require(raw["same_chi_comparison"] is True and raw["stationarity_evaluated"] is True
            and raw["stationarity_passed"] is passed, "worker stationarity decision mismatch")
    before = refinement.interval(rows["fixed_q1"], rows["fixed_q3"], old)
    after = refinement.interval(rows["fixed_q1"], rows["fixed_q3"], new)
    columns = {"fixed_Q1_to_Q3": campaign._signed_column_delta(geometry, rows["fixed_q1"]["sz_profile"],
                                                             rows["fixed_q3"]["sz_profile"], 2),
               "Q3_to_Q5_before": campaign._signed_column_delta(geometry, rows["fixed_q3"]["sz_profile"], old["sz_profile"], 2),
               "Q3_to_Q5_after": campaign._signed_column_delta(geometry, rows["fixed_q3"]["sz_profile"], new["sz_profile"], 2),
               "Q5_10_minus_8": campaign._signed_column_delta(geometry, old["sz_profile"], new["sz_profile"], 0)}
    require(all(item["sum_rule_passed"] for item in columns.values()), "column charge sum rule failed")
    # This shared helper's interpretation text mentions higher/lower Q; the
    # same-Q refinement uses its unchanged signed sum rule with delta_q=0.
    columns["Q5_10_minus_8"]["normalization"] = "column sums of Sz(Q5,cumulative10)-Sz(Q5,cumulative8); no clipping or renormalization"
    aligned = campaign._circumference_shift_comparison(geometry, old, new, 0.0)
    require(aligned["status"] == "computed", "circumference comparison unavailable")
    source_paths = [Path(__file__).resolve(), Path(refinement.__file__).resolve(), Path(campaign.__file__).resolve()]
    return {"schema_version": 1, "recorded_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
            "study": "NN N27 Q5 chi256 cumulative8_to10; Q1/Q3 cumulative8 fixed",
            "claims": {"evidence": "finite variational trials", "rigorous_boundary_error_bounds": False,
                       "matched_sweep_counts": False, "matched_accuracy_inferred": False,
                       "ground_state_convergence": "not_established", "bulk_plateau": "not_established",
                       "phase_identification": "not_attempted", "unexplored_sectors": "Q>=7",
                       "variance": "Hamiltonian dispersion, not a ground-energy error bound"},
            "verification_scope": "saved scalar/profile arithmetic, checkpoint byte hashes and metadata, archived analysis sources; no MPS load/contraction or ED; backend digest maps compared but backend files not rehashed",
            "analysis_source_sha256": {label(path): digest(path) for path in source_paths},
            "inputs": rows, "new_run_verification": verified, "solver": reference["record"]["solver"],
            "runtime": reference["record"]["runtime"], "geometry": geometry, "stationarity_limits": LIMITS,
            "field_convention": "F_Q=E_Q-h*Q/2; h_lower=E3-E1 fixed; h_upper=E5-E3; total exchange energies, J=1",
            "field_intervals": {"before": before, "after": after,
                                "change": {key: after[key] - before[key] for key in ("h_lower", "h_upper", "width")}},
            "refinement": {"Q": 5, "from_cumulative_sweeps": 8, "to_cumulative_sweeps": 10,
                           "energy_change": delta_e, "energy_change_per_site": delta_e / 27,
                           "variance_change": new["variance"] - old["variance"],
                           "last_truncation_change": new["last_sweep_truncation_error"] - old["last_sweep_truncation_error"],
                           "last_sweep_energy_change_recomputed": last_sweep,
                           "scalar_change_convention": "new minus old",
                           "profile_difference_convention": "old minus new at matching coordinates; signed mean follows this convention",
                           "raw_sz_profile_difference": [a - b for a, b in zip(old["sz_profile"], new["sz_profile"])],
                           "raw_bond_profile_difference": [a - b for a, b in zip(old["bond_energy"], new["bond_energy"])],
                           "sz_profile_difference": campaign._difference(old["sz_profile"], new["sz_profile"]),
                           "bond_profile_difference": campaign._difference(old["bond_energy"], new["bond_energy"]),
                           "stationarity_criteria": criteria, "stationarity_passed_recomputed": passed,
                           "y_shift_only_comparison": aligned,
                           "spatial_change_scope": "profiles measured at cumulative8 and cumulative10 only; no attribution to sweep9-to10 alone"},
            "signed_column_changes": columns}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    for role in ROLES:
        parser.add_argument("--" + role.replace("_", "-"), type=Path, required=True)
    parser.add_argument("--new-execution", type=Path, required=True)
    args = parser.parse_args()
    output = args.output.resolve()
    temporary = None
    try:
        require(output.suffix == ".json" and not output.exists(), "use a new .json output path")
        result = analyze({role: getattr(args, role) for role in ROLES}, args.new_execution)
        encoded = json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
        output.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=output.parent, delete=False) as handle:
            temporary = Path(handle.name)
            handle.write(encoded)
        os.link(temporary, output)  # Atomic publication that also refuses a racing overwrite.
    except (OSError, KeyError, TypeError, ValueError) as error:
        parser.exit(2, f"Q5 fixed-chi analysis: {error}\n")
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    print(json.dumps({"output": str(output), "field_intervals": result["field_intervals"]}, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
