#!/usr/bin/env python3
"""Audit one-sweep Q5 chi256/512 forks of the same N27 cumulative8 state.

Finite variational responses and trial-sector intervals only; no MPS contraction,
converged chi extrapolation, rigorous energy error bound, or phase identification.
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
import analyze_q5_fixed_chi as fixed
import analyze_sector_refinement as refinement
import summarize_static_campaign as campaign

require, digest, label, rootpath = refinement.require, refinement.digest, refinement.label, refinement.rootpath
LIMITS = refinement.LIMITS
ROLES = ("fixed_q1", "fixed_q3", "parent_q5", "chi256_q5", "chi512_q5")
RUNTIME_KEYS = ("ITensorMPS", "ITensors", "KrylovKit", "LinearAlgebra", "NDTensors", "SHA",
                "Serialization", "TOML", "arch", "blas_threads", "julia", "julia_threads", "kernel", "word_size")


def load_fork(path: Path, chi: int) -> dict:
    path = path.resolve()
    payload = path.read_bytes()
    record = tomllib.loads(payload.decode("utf-8"))
    name = f"chi{chi}_q5"
    require(record.get("status") == fixed.COMPLETE, f"{name}: worker is not finalized and complete")
    require(record.get("schema_version") == 1 and record.get("N") == 27 and record.get("Q") == 5,
            f"{name}: wrong schema/N/Q")
    require(all(record.get(k) is True for k in ("sources_unchanged", "backend_unchanged", "parent_unchanged")),
            f"{name}: source/parent finalization absent or failed")
    geometry = campaign._geometry(record["configuration"], record)
    require((geometry["Lx"], geometry["Ly"]) == (3, 3) and geometry.get("gauge") == "seam"
            and geometry.get("ordering") == "x_then_y_then_A_B_C", f"{name}: wrong geometry")
    cfg, solver = record["config"], record["solver"]
    require(all(cfg[k] == v for k, v in {"Lx": 3, "Ly": 3, "Q": 5, "maxdim": chi,
            "initialization": "resume", "initial_completed_sweeps": 8, "batches": 1, "batch_sweeps": 1}.items())
            and cfg["stationarity"] == LIMITS, f"{name}: wrong configuration")
    require(solver["maxdim"] == chi and solver["nsweeps"] == 1 and solver["cutoff"] == solver["noise"] == 0,
            f"{name}: wrong solver")
    rows = record["batches"]
    require(len(rows) == 1 and campaign._eligible_batch(rows[0]) and record["completed_sweeps"] == 9,
            f"{name}: no unique completed cumulative9 batch")
    row = rows[0]
    require(row["cumulative_sweeps"] == 9 and row["batch"] == 1 and row["maxlinkdim"] == chi,
            f"{name}: wrong sweep count or actual chi")
    require(row.get("within_batch_sweep_change_measured") is False and "last_sweep_energy_change" not in row,
            f"{name}: one sweep cannot supply a final two-sweep difference")
    energy, variance = (campaign._number(row[k], k) for k in ("energy", "variance"))
    sz = campaign._vector(row["sz_profile"], 27, "sz_profile")
    bonds = campaign._vector(row["bond_energy"], len(geometry["bonds"]), "bond_energy")
    truncation = campaign._vector(row["measured_truncation_errors"], 1, "truncation")
    energies = campaign._vector(row["sweep_energies"], 1, "sweep_energies")
    require(0 <= truncation[0] <= 1 and abs(math.fsum(sz) - 2.5) <= 1e-10
            and abs(math.fsum(bonds) - energy) <= 1e-8, f"{name}: observable sum rule failed")
    require(variance >= -100 * sys.float_info.epsilon * max(1, energy * energy),
            f"{name}: negative variance beyond roundoff")
    cp = rootpath(row["checkpoint"])
    hashes = {name: digest(cp / name) for name in ("metadata.toml", "state.jls", "checksums.toml")}
    require(hashes["metadata.toml"] == row["checkpoint_metadata_sha256"]
            and hashes["state.jls"] == row["checkpoint_payload_sha256"], f"{name}: checkpoint digest mismatch")
    checksums = tomllib.loads((cp / "checksums.toml").read_text())
    for filename in ("metadata.toml", "state.jls"):
        require(checksums["files"][filename] == {"sha256": hashes[filename], "bytes": (cp / filename).stat().st_size},
                f"{name}: checkpoint checksum table mismatch")
    meta = tomllib.loads((cp / "metadata.toml").read_text())
    require(meta["status"] == "trial" and meta["theta_path"] == [0.0]
            and meta["configuration"] == record["configuration"], f"{name}: checkpoint configuration mismatch")
    state = meta["state"]
    require(state["Q"] == 5 and state["theta"] == 0 and state["energy"] == energy and state["sz"] == sz
            and state["sweep_energies"] == energies and state["max_truncation_errors"] == truncation,
            f"{name}: checkpoint measurements mismatch")
    backend = record["code"]["source_sha256"]
    require(isinstance(backend, dict) and bool(backend) and meta["provenance"]["source_sha256"] == backend,
            f"{name}: checkpoint/backend identity mismatch")
    require(meta["runtime"] == {k: v for k, v in record["runtime"].items()
                                if k not in ("julia_threads", "blas_threads")}, f"{name}: checkpoint runtime mismatch")
    summary = {"validation_path": label(path), "validation_sha256": hashlib.sha256(payload).hexdigest(),
               "case_id": record["case_id"], "N": 27, "Q": 5, "cumulative_sweeps": 9,
               "configured_maxdim": chi, "actual_maxlinkdim": chi, "selected_batch": 1,
               "energy": energy, "variance": variance, "variance_per_site": variance / 27,
               "last_sweep_truncation_error": truncation[0], "sweep_energies": energies,
               "sz_profile": sz, "bond_energy": bonds, "checkpoint": label(cp), "checkpoint_sha256": hashes,
               "source_finalization": {k: record[k] for k in ("sources_unchanged", "backend_unchanged", "parent_unchanged")},
               "backend_source_sha256": backend, "recorded_analysis_source_sha256": record["analysis_source_sha256"],
               "quality": campaign._quality_summary(row)}
    return {"path": path, "record": record, "row": row, "geometry": geometry, "summary": summary}


def check_execution_and_archive(case: dict, *, historical_q3: bool = False) -> dict:
    record, summary = case["record"], case["summary"]
    execution_path = campaign._companion_paths(case["path"])[0]
    execution = tomllib.loads(execution_path.read_text())
    require(execution.get("status") == "exited" and execution.get("worker_exit_confirmed") is True
            and execution.get("worker_exit_code") == 0, f"{summary['case_id']}: worker exit unconfirmed")
    require(0 < execution["wall_limit_seconds"] <= 600
            and 0 <= execution["wall_elapsed_seconds"] <= execution["wall_limit_seconds"]
            and record["wall_limit_seconds"] == execution["wall_limit_seconds"], "worker exceeded its wall budget")
    require(execution["julia_threads"] == execution["blas_threads"] == 1
            and record["runtime"]["julia_threads"] == record["runtime"]["blas_threads"] == 1, "wrong thread allocation")
    run_root = rootpath(summary["checkpoint"]).parent.parent
    require(digest(run_root / "validation.toml") == summary["validation_sha256"]
            and digest(run_root / "execution.toml") == digest(execution_path), "run/archive binding mismatch")
    settings = tomllib.loads((rootpath(summary["checkpoint"]) / "metadata.toml").read_text())["settings"]
    for key, value in record["solver"].items():
        require(settings[key] == ([value] if key == "maxdim" else value), f"checkpoint solver mismatch: {key}")
    archived = run_root / "analysis-sources"
    result = {"execution_path": label(execution_path), "execution_sha256": digest(execution_path),
              "execution": campaign._execution_summary(execution), "checkpoint_solver_verified": True}
    config_archive = archived / "config.toml"
    if historical_q3 and not config_archive.exists():
        config_archive = archived / "examples/configs/nn27_q3_chi256_resume.toml"
    require(digest(config_archive) == record["config_sha256"]
            and tomllib.loads(config_archive.read_text()) == record["config"], "archived configuration mismatch")
    require(record["analysis_source_sha256"] and all(digest(archived / name) == sha
            for name, sha in record["analysis_source_sha256"].items()), "archived analysis sources mismatch")
    return {**result, "archived_configuration_and_sources_verified": True}


def check_parent(parent_case: dict, fork: dict) -> dict:
    before, record = parent_case["summary"], fork["record"]
    cfg, parent = record["config"], record["parent"]
    require(parent.get("strict_source_runtime_configuration_load") is True
            and parent["initial_completed_sweeps"] == 8, "strict cumulative8 parent load absent")
    cp, evidence = rootpath(before["checkpoint"]), rootpath(parent["record"])
    require(rootpath(cfg["parent_checkpoint"]) == rootpath(parent["checkpoint"]) == cp, "parent checkpoint mismatch")
    require(rootpath(cfg["parent_record"]) == evidence and digest(evidence) == before["validation_sha256"]
            == parent["record_sha256"] == cfg["parent_record_sha256"], "parent evidence digest mismatch")
    require(cfg["parent_metadata_sha256"] == before["checkpoint_sha256"]["metadata.toml"], "parent metadata mismatch")
    expected = {cp / name: sha for name, sha in before["checkpoint_sha256"].items()}
    expected[evidence] = before["validation_sha256"]
    require({rootpath(k): v for k, v in parent["checkpoint_sha256"].items()} == expected, "parent digest table mismatch")
    require(abs(parent["loaded_energy"] - before["energy"]) <= 1e-12, "loaded parent energy mismatch")
    return {"status": "verified_same_cumulative8_parent", "checkpoint": label(cp),
            "parent_validation_sha256": before["validation_sha256"]}


def compare(geometry: dict, left: dict, right: dict) -> dict:
    aligned = campaign._circumference_shift_comparison(geometry, left, right, 0.0)
    require(aligned["status"] == "computed", "circumference comparison unavailable")
    column = campaign._signed_column_delta(geometry, left["sz_profile"], right["sz_profile"], 0)
    require(column["sum_rule_passed"], "same-sector column sum rule failed")
    column["normalization"] = "column sums of right Sz minus left Sz; no clipping or renormalization"
    return {"scalar_change_convention": "right minus left", "energy_change": right["energy"] - left["energy"],
            "energy_change_per_site": (right["energy"] - left["energy"]) / 27,
            "variance_change": right["variance"] - left["variance"],
            "last_truncation_change": right["last_sweep_truncation_error"] - left["last_sweep_truncation_error"],
            "profile_difference_convention": "left minus right; signed mean follows this convention",
            "sz_profile_difference": campaign._difference(left["sz_profile"], right["sz_profile"]),
            "bond_profile_difference": campaign._difference(left["bond_energy"], right["bond_energy"]),
            "y_shift_only_comparison": aligned, "signed_column_change": column}


def parent_criteria(parent: dict, fork: dict, comparison: dict) -> dict:
    row, result = fork["row"], fork["summary"]
    same_chi = result["configured_maxdim"] == parent["configured_maxdim"]
    optimizer_error = abs(result["energy"] - result["sweep_energies"][-1])
    values = {"energy_per_site": max(abs(comparison["energy_change"]), optimizer_error) / 27,
              "sz_profile": max(abs(a - b) for a, b in zip(parent["sz_profile"], result["sz_profile"])),
              "bond_profile": max(abs(a - b) for a, b in zip(parent["bond_energy"], result["bond_energy"])),
              "variance_per_site": abs(result["variance"]) / 27, "truncation": result["last_sweep_truncation_error"]}
    criteria = {k: {"value": v, "limit": LIMITS[k], "passed": v <= LIMITS[k]} for k, v in values.items()}
    passed = all(item["passed"] for item in criteria.values())
    for key, value in {"energy_change_from_previous": comparison["energy_change"],
                       "last_optimizer_energy_error": optimizer_error, "max_sz_change": values["sz_profile"],
                       "max_bond_change": values["bond_profile"], "variance_per_site": result["variance"] / 27}.items():
        require(abs(campaign._number(row[key], key) - value) <= 5e-13, f"worker arithmetic mismatch: {key}")
    require(row["same_chi_comparison"] is same_chi and row["stationarity_evaluated"] is same_chi
            and row["stationarity_passed"] is (same_chi and passed)
            and row["comparison_kind"] == ("same_chi" if same_chi else "chi_change"), "worker stationarity mismatch")
    return {"comparison_kind": "same_chi_one_sweep" if same_chi else "chi_change_one_sweep",
            "individual_precision_thresholds": criteria, "last_optimizer_energy_error": optimizer_error,
            "within_batch_sweep_change_measured": False, "final_two_sweep_difference": None,
            "stationarity_evaluated": same_chi, "stationarity_passed": passed if same_chi else None,
            "interpretation": "one parent-to-fork response; no final two-sweep difference; chi changes are not stationarity tests"}


def analyze(paths: dict[str, Path]) -> dict:
    cases = {role: load_fork(paths[role], chi) for role, chi in (("chi256_q5", 256), ("chi512_q5", 512))}
    for role, fixed_role in (("fixed_q1", "fixed_q1"), ("fixed_q3", "fixed_q3"), ("parent_q5", "old_q5")):
        cases[role] = fixed.load_input(paths[role], fixed_role)
    if "previous25610" in paths:
        cases["previous25610"] = fixed.load_input(paths["previous25610"], "new_q5")
    reference = cases["parent_q5"]
    geometry = {k: v for k, v in reference["geometry"].items() if k != "Q"}
    for role, case in cases.items():
        require({k: v for k, v in case["geometry"].items() if k != "Q"} == geometry, f"{role}: geometry mismatch")
        require(case["record"]["runtime"] == reference["record"]["runtime"], f"{role}: runtime mismatch")
        require(case["summary"]["backend_source_sha256"] == reference["summary"]["backend_source_sha256"],
                f"{role}: backend source mismatch")
        require({k: v for k, v in case["record"]["solver"].items() if k not in ("maxdim", "nsweeps")}
                == {k: v for k, v in reference["record"]["solver"].items() if k not in ("maxdim", "nsweeps")},
                f"{role}: solver mismatch beyond maxdim and historical batch length")
    a, b = (cases[key]["record"] for key in ("chi256_q5", "chi512_q5"))
    require({k: v for k, v in a["solver"].items() if k != "maxdim"}
            == {k: v for k, v in b["solver"].items() if k != "maxdim"}, "fork solvers differ beyond maxdim")
    require({k: v for k, v in a["config"].items() if k not in ("maxdim", "case_id")}
            == {k: v for k, v in b["config"].items() if k not in ("maxdim", "case_id")}, "fork configs differ beyond maxdim/case_id")
    require(a["analysis_source_sha256"] == b["analysis_source_sha256"], "fork analysis source mismatch")
    verifications = {role: check_execution_and_archive(case, historical_q3=role == "fixed_q3")
                     for role, case in cases.items()}
    rows = {role: case["summary"] for role, case in cases.items()}
    comparisons = {}
    for role in ("chi256_q5", "chi512_q5"):
        verifications[role]["parent_link"] = check_parent(reference, cases[role])
        comparison = compare(geometry, rows["parent_q5"], rows[role])
        comparison["precision"] = parent_criteria(rows["parent_q5"], cases[role], comparison)
        comparisons[f"parent_to_{role}"] = comparison
    comparisons["chi256_to_chi512"] = compare(geometry, rows["chi256_q5"], rows["chi512_q5"])
    comparisons["chi256_to_chi512"]["interpretation"] = "controlled one-sweep solver response from a common parent; not converged chi extrapolation or stationarity"
    q5_roles = ["parent_q5", "chi256_q5", "chi512_q5"]
    if "previous25610" in cases:
        q5_roles.append("previous25610")
        previous = cases["previous25610"]
        verifications["previous25610"]["parent_link"] = fixed.check_new_run(
            reference, previous, campaign._companion_paths(previous["path"])[0])
        comparison = compare(geometry, rows["chi256_q5"], rows["previous25610"])
        comparison["interpretation"] = "same chi with one additional nominal sweep in a separate run from the verified common cumulative8 parent; descriptive difference, not a direct continuation of the new cumulative9 checkpoint"
        expectation = rows["chi256_q5"]["energy"]
        optimizer = rows["previous25610"]["sweep_energies"][0]
        comparison["cumulative9_energy_agreement_diagnostic"] = {
            "one_sweep_run_final_MPO_expectation_energy": expectation,
            "two_sweep_run_first_reported_solver_energy": optimizer,
            "MPO_expectation_minus_reported_solver_energy": expectation - optimizer,
            "absolute_difference": abs(expectation - optimizer),
            "interpretation": "energy agreement diagnostic across separate runs and measurement types; no bitwise equality required or wavefunction identity inferred",
            "new_one_sweep_run_final_two_sweep_difference": None}
        comparisons["chi256_9_to_previous256_10"] = comparison
        comparisons["chi512_9_to_previous256_10"] = compare(geometry, rows["chi512_q5"], rows["previous25610"])
        comparisons["chi512_9_to_previous256_10"]["interpretation"] = "unmatched chi and sweep count; descriptive comparison only"
    intervals = {role: refinement.interval(rows["fixed_q1"], rows["fixed_q3"], rows[role]) for role in q5_roles}
    columns = {"fixed_Q1_to_Q3": campaign._signed_column_delta(geometry, rows["fixed_q1"]["sz_profile"],
                                                              rows["fixed_q3"]["sz_profile"], 2)}
    columns.update({f"Q3_to_{role}": campaign._signed_column_delta(geometry, rows["fixed_q3"]["sz_profile"],
                                                                  rows[role]["sz_profile"], 2) for role in q5_roles})
    require(all(item["sum_rule_passed"] for item in columns.values()), "sector column charge sum rule failed")
    return {"schema_version": 1, "recorded_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
            "study": "NN N27 Q5 matched chi256/512 one-sweep cumulative8-to9 forks; Q1/Q3 cumulative8 fixed",
            "claims": {"evidence": "finite variational trials", "matched_q5_parent_and_sweep_count": True,
                       "matched_accuracy_inferred": False, "rigorous_boundary_error_bounds": False,
                       "converged_chi_extrapolation": False, "ground_state_convergence": "not_established",
                       "bulk_plateau": "not_established", "phase_identification": "not_attempted",
                       "unexplored_sectors": "Q>=7; sector skipping not excluded",
                       "variance": "Hamiltonian dispersion, not a ground-energy error bound"},
            "verification_scope": "saved measurements, checkpoint hashes and metadata, execution, archived config/drivers; no MPS load/contraction or ED; backend maps compared without rehashing backend files; historical Q3 backend_unchanged absence retained",
            "analysis_source_sha256": {label(Path(module.__file__).resolve()): digest(Path(module.__file__).resolve())
                                       for module in (sys.modules[__name__], fixed, refinement, campaign)},
            "inputs": rows, "run_verification": verifications, "geometry": geometry,
            "runtime": {k: reference["record"]["runtime"][k] for k in RUNTIME_KEYS if k in reference["record"]["runtime"]},
            "fork_solvers": {role: cases[role]["record"]["solver"] for role in ("chi256_q5", "chi512_q5")},
            "stationarity_limits": LIMITS, "comparisons": comparisons,
            "field_convention": "F_Q=E_Q-h*Q/2; h_lower=E3-E1 fixed; h_upper=E5-E3; total exchange energies, J=1",
            "field_intervals": intervals, "signed_sector_column_changes": columns}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    for role in ROLES:
        parser.add_argument("--" + role.replace("_", "-"), type=Path, required=True)
    parser.add_argument("--previous25610", type=Path)
    args = parser.parse_args()
    output, temporary = args.output.resolve(), None
    try:
        require(output.suffix == ".json" and not output.exists(), "use a new .json output path")
        paths = {role: getattr(args, role) for role in ROLES}
        if args.previous25610 is not None:
            paths["previous25610"] = args.previous25610
        result = analyze(paths)
        encoded = json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
        output.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=output.parent, delete=False) as handle:
            temporary = Path(handle.name)
            handle.write(encoded)
        os.link(temporary, output)
    except (OSError, KeyError, TypeError, ValueError) as error:
        parser.exit(2, f"Q5 matched-chi analysis: {error}\n")
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    print(json.dumps({"output": str(output), "field_intervals": result["field_intervals"]}, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
