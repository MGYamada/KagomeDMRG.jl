#!/usr/bin/env python3
"""Audit bounded NN N27 chi512 follow-ups using saved trials only.

Q1 changes chi and sweep count together. Q5 alone provides a same-chi
parent-to-child stationarity check. Neither comparison establishes matched
accuracy or a converged field boundary; no MPS is contracted here.
"""
from __future__ import annotations

import argparse
import datetime as dt
import json
import math
import os
from pathlib import Path
import sys
import tempfile
import tomllib

sys.dont_write_bytecode = True
import analyze_q5_matched_chi as matched
import analyze_sector_refinement as refinement
import summarize_static_campaign as campaign

require, digest, label, rootpath = refinement.require, refinement.digest, refinement.label, refinement.rootpath
LIMITS = refinement.LIMITS
DATA = Path(__file__).resolve().parents[1] / "docs/research/data"


def toml(path: Path) -> dict:
    return tomllib.loads(path.read_text())


def checkpoint_summary(path: Path, record: dict, row: dict, q: int, chi: int, sweeps: int) -> dict:
    """Validate a completed row and its immutable checkpoint independently of worker completion."""
    geometry = campaign._geometry(record["configuration"], record)
    require(record["N"] == 27 and record["Q"] == q and (geometry["Lx"], geometry["Ly"]) == (3, 3)
            and geometry["gauge"] == "seam" and geometry["ordering"] == "x_then_y_then_A_B_C", "wrong geometry/Q")
    require(campaign._eligible_batch(row) and row["cumulative_sweeps"] == sweeps
            and 1 <= row["maxlinkdim"] <= chi and record["solver"]["maxdim"] == chi, "ineligible selected diagnostic row")
    energy, variance = (campaign._number(row[k], k) for k in ("energy", "variance"))
    nsweeps = record["solver"]["nsweeps"]
    sz = campaign._vector(row["sz_profile"], 27, "Sz")
    bonds = campaign._vector(row["bond_energy"], len(geometry["bonds"]), "bonds")
    trunc = campaign._vector(row["measured_truncation_errors"], nsweeps, "truncation")
    energies = campaign._vector(row["sweep_energies"], nsweeps, "sweep energies")
    require(abs(math.fsum(sz) - q / 2) <= 1e-10 and abs(math.fsum(bonds) - energy) <= 1e-8
            and all(0 <= x <= 1 for x in trunc), "observable sum rule/truncation failed")
    scale = 100 * sys.float_info.epsilon * max(1, energy * energy)
    require(variance >= -scale, "negative variance beyond roundoff")
    moments = {"raw_complex_second_moment_recorded": all(k in row for k in ("HdaggerH_real", "HdaggerH_imag"))}
    if moments["raw_complex_second_moment_recorded"]:
        real, imag = (campaign._number(row[k], k) for k in ("HdaggerH_real", "HdaggerH_imag"))
        require(abs(imag) <= scale and abs(real - energy * energy - variance) <= scale, "second moment arithmetic failed")
        moments.update(HdaggerH_real=real, HdaggerH_imag=imag, variance_arithmetic_passed=True)
    cp = rootpath(row["checkpoint"])
    hashes = {name: digest(cp / name) for name in ("metadata.toml", "state.jls", "checksums.toml")}
    require(hashes["metadata.toml"] == row["checkpoint_metadata_sha256"]
            and hashes["state.jls"] == row["checkpoint_payload_sha256"], "checkpoint digest mismatch")
    checksums, meta = toml(cp / "checksums.toml"), toml(cp / "metadata.toml")
    for name in ("metadata.toml", "state.jls"):
        require(checksums["files"][name] == {"sha256": hashes[name], "bytes": (cp / name).stat().st_size}, "checkpoint checksum table mismatch")
    state = meta["state"]
    require(meta["status"] == "trial" and meta["theta_path"] == [0.0]
            and meta["configuration"] == record["configuration"], "checkpoint configuration mismatch")
    require(state["Q"] == q and state["theta"] == 0 and state["energy"] == energy and state["sz"] == sz
            and state["sweep_energies"] == energies and state["max_truncation_errors"] == trunc, "checkpoint measurements mismatch")
    backend = record["code"]["source_sha256"]
    require(backend and meta["provenance"]["source_sha256"] == backend, "checkpoint backend mismatch")
    require(meta["runtime"] == {k: v for k, v in record["runtime"].items() if k not in ("julia_threads", "blas_threads")}, "checkpoint runtime mismatch")
    for key, value in record["solver"].items():
        require(meta["settings"][key] == ([value] if key == "maxdim" else value), f"checkpoint solver mismatch: {key}")
    return {"path": path, "record": record, "row": row, "geometry": geometry,
            "summary": {"validation_path": label(path), "validation_sha256": digest(path),
                        "case_id": record["case_id"], "N": 27, "Q": q, "cumulative_sweeps": sweeps,
                        "configured_maxdim": chi, "actual_maxlinkdim": row["maxlinkdim"],
                        "selected_batch": row["batch"], "energy": energy, "variance": variance,
                        "variance_per_site": variance / 27, "last_sweep_truncation_error": trunc[-1],
                        "sz_profile": sz, "bond_energy": bonds, "sweep_energies": energies,
                        "checkpoint": label(cp), "checkpoint_sha256": hashes,
                        "backend_source_sha256": backend, "quality": campaign._quality_summary(row),
                        "moment_verification": moments, "worker_status": record["status"],
                        "source_finalization": {k: record.get(k) for k in ("sources_unchanged", "backend_unchanged", "parent_unchanged")}}}


def archive_check(path: Path, record: dict, run_root: Path, *, limit: int, timed_out: bool = False) -> dict:
    execution_path = campaign._companion_paths(path)[0]
    execution = toml(execution_path)
    require(execution["worker_exit_confirmed"] is True, "worker exit unconfirmed")
    require(execution["status"] == ("timed_out" if timed_out else "exited"), "unexpected execution status")
    require((execution["worker_exit_code"] != 0 if timed_out else execution["worker_exit_code"] == 0), "unexpected worker exit code")
    require(0 < execution["wall_limit_seconds"] <= limit and 0 <= execution["wall_elapsed_seconds"] <= execution["wall_limit_seconds"]
            and record["wall_limit_seconds"] == execution["wall_limit_seconds"], "worker budget mismatch")
    require(execution["julia_threads"] == execution["blas_threads"] == 1
            and record["runtime"]["julia_threads"] == record["runtime"]["blas_threads"] == 1, "thread allocation mismatch")
    require(digest(run_root / "validation.toml") == digest(path)
            and digest(run_root / "execution.toml") == digest(execution_path), "archive/run mismatch")
    archive = run_root / "analysis-sources"
    config = archive / "config.toml"
    if not config.exists():
        config = archive / record["config_path"]
    require(digest(config) == record["config_sha256"] and toml(config) == record["config"], "archived configuration mismatch")
    require(bool(record["analysis_source_sha256"]) and all(digest(archive / p) == h for p, h in record["analysis_source_sha256"].items()), "archived source mismatch")
    return {"execution_path": label(execution_path), "execution_sha256": digest(execution_path),
            "execution": campaign._execution_summary(execution), "archived_config_and_sources_verified": True}


def historical(path: Path, q: int, chi: int, sweeps: int, *, partial_audit: Path | None = None) -> dict:
    path = path.resolve()
    record = toml(path)
    if partial_audit is None:
        require(record["status"] == matched.fixed.COMPLETE and record.get("sources_unchanged") is True
                and record.get("parent_unchanged") is True, "historical worker not finalized")
        require(record.get("backend_unchanged") is True or (q == 3 and chi == 256 and "backend_unchanged" not in record), "historical backend finalization absent")
    else:
        audit = json.loads(partial_audit.read_text())
        require(record["status"] == "running" and audit["status"] == "passed" and all(v is True for v in audit["checks"].values())
                and audit["completed_diagnostic_sweeps"] == [sweeps] and audit["worker_finalization"] == "interrupted_during_second_sweep"
                and audit["validation_sha256"] == digest(path), "historical interrupted worker audit mismatch")
    chosen = [r for r in record["batches"] if campaign._eligible_batch(r) and r["cumulative_sweeps"] == sweeps]
    require(len(chosen) == 1, "required historical completed row absent")
    case = checkpoint_summary(path, record, chosen[0], q, chi, sweeps)
    verify = archive_check(path, record, rootpath(case["summary"]["checkpoint"]).parent.parent,
                           limit=600, timed_out=partial_audit is not None)
    if partial_audit is not None:
        require(audit["execution_sha256"] == verify["execution_sha256"], "post-run audit execution mismatch")
        verify["post_run_hash_audit"] = {"path": label(partial_audit), "sha256": digest(partial_audit),
                                         "selected_completed_sweeps": sweeps, "whole_worker_completed": False,
                                         "scope": "completed first batch only; second batch interrupted; saved audit retained without retroactive worker finalization"}
    case["verification"] = verify
    return case


def parent_check(parent: dict, child: dict) -> dict:
    before, record = parent["summary"], child["record"]
    cfg, link = record["config"], record["parent"]
    require(link["strict_source_runtime_configuration_load"] is True and link["initial_completed_sweeps"]
            == cfg["initial_completed_sweeps"] == before["cumulative_sweeps"], "strict parent sweep link mismatch")
    cp, evidence = rootpath(before["checkpoint"]), rootpath(link["record"])
    require(rootpath(cfg["parent_checkpoint"]) == rootpath(link["checkpoint"]) == cp
            and rootpath(cfg["parent_record"]) == evidence, "parent paths mismatch")
    require(digest(evidence) == link["record_sha256"] == cfg["parent_record_sha256"] == before["validation_sha256"]
            and cfg["parent_metadata_sha256"] == before["checkpoint_sha256"]["metadata.toml"], "parent digest mismatch")
    expected = {cp / k: v for k, v in before["checkpoint_sha256"].items()}
    expected[evidence] = before["validation_sha256"]
    require({rootpath(k): v for k, v in link["checkpoint_sha256"].items()} == expected
            and abs(link["loaded_energy"] - before["energy"]) <= 1e-12, "parent identity mismatch")
    return {"strict_parent_link_verified": True, "parent_validation_sha256": before["validation_sha256"],
            "parent_cumulative_sweeps": before["cumulative_sweeps"]}


def split_case(path: Path, q: int, sweeps: int) -> dict:
    path = path.resolve()
    record = toml(path)
    require(record.get("stage") == "diagnose" and record.get("status") == "completed_saved_state_diagnostics_accuracy_separate",
            "a completed diagnostic record is required; solve-only or interrupted records are insufficient")
    flags = ("sources_unchanged", "backend_unchanged", "parent_unchanged", "inputs_unchanged")
    require(all(record.get(k) is True for k in flags), "diagnostic provenance finalization failed")
    cfg = record["config"]
    require(all(cfg[k] == v for k, v in {"Lx": 3, "Ly": 3, "Q": q, "maxdim": 512,
            "initialization": "resume", "initial_completed_sweeps": sweeps - 1, "batches": 1, "batch_sweeps": 1}.items())
            and cfg["stationarity"] == LIMITS, "wrong split-run configuration")
    require(record["completed_sweeps"] == sweeps and len(record["batches"]) == 1, "wrong diagnostic sweep count")
    row = record["batches"][0]
    require(row["saved_state_remeasurement_passed"] is True and row["within_batch_sweep_change_measured"] is False
            and "last_sweep_energy_change" not in row, "remeasurement missing or invalid single-sweep claim")
    require(0 <= campaign._number(row["remeasured_energy_error"], "energy error") <= 1e-10 * max(1, abs(row["energy"]))
            and all(0 <= campaign._number(row[k], k) <= 1e-12 for k in ("remeasured_sz_error", "remeasured_norm_error", "remeasured_charge_error")),
            "saved-state remeasurement exceeds tolerance")
    case = checkpoint_summary(path, record, row, q, 512, sweeps)
    solve_path, execution_path = rootpath(record["solve_record"]), rootpath(record["solve_execution"])
    require(digest(solve_path) == record["solve_record_sha256"] and digest(execution_path) == record["solve_execution_sha256"], "diagnostic solve binding mismatch")
    solve = toml(solve_path)
    require(solve["stage"] == "solve" and solve["status"] == "completed_solve_diagnostics_deferred"
            and all(solve.get(k) is True for k in flags), "solve has not finalized successfully")
    for key in ("config", "config_sha256", "configuration", "solver", "runtime", "parent", "analysis_source_sha256"):
        require(solve[key] == record[key], f"solve/diagnostic mismatch: {key}")
    require(solve["code"]["source_sha256"] == record["code"]["source_sha256"], "solve/diagnostic backend mismatch")
    require(len(solve["batches"]) == 1 and solve["completed_sweeps"] == sweeps, "wrong solve row count")
    original = solve["batches"][0]
    require(original["status"] == "checkpoint_saved_diagnostics_deferred" and original["checkpoint_verified"] is True
            and all(k not in original for k in ("integrity_passed", "variance", "last_sweep_energy_change", "stationarity_passed")),
            "solve-only record claims unperformed diagnostics")
    for key in ("checkpoint", "checkpoint_sha256", "checkpoint_metadata_sha256", "checkpoint_payload_sha256", "settings",
                "energy", "sz_profile", "sweep_energies", "measured_truncation_errors", "maxlinkdim", "cumulative_sweeps"):
        require(row[key] == original[key], f"diagnostics altered solve measurement: {key}")
    cp = rootpath(row["checkpoint"])
    require(cp.parent.parent == solve_path.parent and execution_path == solve_path.parent / "execution.toml", "solve checkpoint/execution ownership mismatch")
    expected_hashes = {cp / k: v for k, v in case["summary"]["checkpoint_sha256"].items()}
    require({rootpath(k): v for k, v in row["checkpoint_sha256"].items()} == expected_hashes
            and toml(cp / "metadata.toml")["settings"] == row["settings"], "solve checkpoint pin/settings mismatch")
    for stage in (solve, record):
        require(stage["input_sha256"] and all(digest(rootpath(p)) == sha for p, sha in stage["input_sha256"].items()), "split-run input bytes changed")
    # The diagnostic manifest may be copied into docs; its original output is
    # identified by the supervised command instead of the solve-owned checkpoint.
    diagnostic_execution = toml(campaign._companion_paths(path)[0])
    command = diagnostic_execution["command"]
    require(len(command) == 11 and Path(command[0]).name == "julia"
            and command[1:4] == ["--project=.", "--startup-file=no", "--threads=1"]
            and rootpath(command[4]) == rootpath("examples/research_static_split.jl") and command[5] == "diagnose"
            and rootpath(command[7]) == rootpath(record["config_path"])
            and rootpath(command[8]) == solve_path.parent and command[9] == digest(solve_path)
            and command[10] == digest(execution_path), "diagnostic command does not bind solve/configuration")
    solve_command = toml(execution_path)["command"]
    require(len(solve_command) == 8 and solve_command[:5] == command[:5] and solve_command[5] == "solve"
            and rootpath(solve_command[6]) == solve_path.parent and solve_command[7] == command[7], "solve command mismatch")
    verification = {"solve": archive_check(solve_path, solve, solve_path.parent, limit=900),
                    "diagnostics": archive_check(path, record, rootpath(command[6]), limit=180)}
    case["verification"] = verification
    verification.update(solve_record_path=label(solve_path), solve_record_sha256=digest(solve_path),
                        solve_only_status_preserved=True, diagnostic_remeasurement_verified=True,
                        immutable_input_hashes_verified=True)
    case["summary"]["parent_overlap_abs"] = campaign._number(row["parent_overlap_abs"], "parent overlap")
    return case


def analyze(q1_path: Path, q5_path: Path) -> dict:
    cases = {
        "q1_256_8": historical(DATA / "p4_sector_refine_nn27_q1_chi256_validation.toml", 1, 256, 8),
        "q3_256_8": historical(DATA / "p4_campaign_nn27_q3_chi256_validation.toml", 3, 256, 8),
        "q5_256_8": historical(DATA / "p4_sector_refine_nn27_q5_chi256_validation.toml", 5, 256, 8),
        "q3_512_9": historical(DATA / "p4_campaign_nn27_q3_chi512_validation.toml", 3, 512, 9,
                                 partial_audit=DATA / "p4_campaign_nn27_q3_chi512_post_run_hash_audit.json"),
        "q5_512_9": historical(DATA / "p4_q5_matched_chi512_8to9_validation.toml", 5, 512, 9),
        "q1_512_9": split_case(q1_path, 1, 9),
        "q5_512_10": split_case(q5_path, 5, 10),
    }
    reference = cases["q1_256_8"]
    geometry = {k: v for k, v in reference["geometry"].items() if k != "Q"}
    rows = {role: case["summary"] for role, case in cases.items()}
    for role, case in cases.items():
        require({k: v for k, v in case["geometry"].items() if k != "Q"} == geometry, f"{role}: geometry mismatch")
        require(case["record"]["runtime"] == reference["record"]["runtime"], f"{role}: runtime mismatch")
        require(rows[role]["backend_source_sha256"] == reference["summary"]["backend_source_sha256"], f"{role}: backend mismatch")
        require({k: v for k, v in case["record"]["solver"].items() if k not in ("nsweeps", "maxdim")}
                == {k: v for k, v in reference["record"]["solver"].items() if k not in ("nsweeps", "maxdim")}, f"{role}: solver mismatch")
    comparisons = {}
    for old, new in (("q3_256_8", "q3_512_9"), ("q5_256_8", "q5_512_9")):
        cases[new]["verification"]["parent_link"] = parent_check(cases[old], cases[new])
    for old, new in (("q1_256_8", "q1_512_9"), ("q5_512_9", "q5_512_10")):
        cases[new]["verification"]["parent_link"] = parent_check(cases[old], cases[new])
        comparison = matched.compare(geometry, rows[old], rows[new])
        comparison["precision"] = matched.parent_criteria(rows[old], cases[new], comparison)
        thresholds = comparison["precision"]["individual_precision_thresholds"]
        raw = cases[new]["row"]
        require(set(raw["precision_values"]) == set(raw["precision_passed"]) == set(LIMITS)
                and all(abs(raw["precision_values"][k] - v["value"]) <= 5e-13
                        and raw["precision_passed"][k] is v["passed"] for k, v in thresholds.items())
                and raw["all_precision_conditions_passed"] is all(v["passed"] for v in thresholds.values()),
                "worker five-threshold arithmetic mismatch")
        comparison["interpretation"] = ("chi and one additional sweep change together; pure chi effect not isolated"
                                          if old.startswith("q1") else "direct same-chi one-sweep continuation; matched accuracy not inferred")
        comparison["from"] = old
        comparison["to"] = new
        comparisons[new] = comparison
    selections = {
        "all256_cumulative8": ("q1_256_8", "q3_256_8", "q5_256_8"),
        "previous_mixed_Q5_512_9": ("q1_256_8", "q3_256_8", "q5_512_9"),
        "replace_Q3_only": ("q1_256_8", "q3_512_9", "q5_512_9"),
        "all512_cumulative9": ("q1_512_9", "q3_512_9", "q5_512_9"),
        "current512_Q5_cumulative10": ("q1_512_9", "q3_512_9", "q5_512_10"),
    }
    intervals, columns, previous = {}, {}, None
    for name, selection in selections.items():
        interval = refinement.interval(*(rows[role] for role in selection))
        interval["selection"] = dict(zip(("Q1", "Q3", "Q5"), selection))
        if previous is not None:
            interval["change_from_previous_selection"] = {key: interval[key] - previous[key] for key in ("h_lower", "h_upper", "width")}
        intervals[name], previous = interval, interval
        columns[name] = {"Q1_to_Q3": campaign._signed_column_delta(geometry, rows[selection[0]]["sz_profile"], rows[selection[1]]["sz_profile"], 2),
                         "Q3_to_Q5": campaign._signed_column_delta(geometry, rows[selection[1]]["sz_profile"], rows[selection[2]]["sz_profile"], 2)}
    require(all(v["sum_rule_passed"] for cols in columns.values() for v in cols.values()), "sector charge sum rule failed")
    modules = (sys.modules[__name__], matched, matched.fixed, refinement, campaign)
    return {"schema_version": 1, "recorded_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
            "study": "NN N27 Q1 chi256 cumulative8 to chi512 cumulative9; Q5 chi512 cumulative9 to10",
            "claims": {"evidence": "finite variational trials", "matched_accuracy_inferred": False,
                       "Q1_pure_chi_effect_isolated": False, "Q5_same_chi_stationarity_evaluated": True,
                       "ground_state_convergence": "not_established", "bulk_plateau": "not_established",
                       "phase_identification": "not_attempted", "rigorous_boundary_error_bounds": False,
                       "unexplored_sectors": "Q>=7; sector skipping not excluded",
                       "variance": "Hamiltonian dispersion; not a rigorous ground-energy error bound"},
            "verification_scope": "saved measurements, strict parent and split-run links, checkpoint byte hashes and metadata, archived config/code, runtime/settings/execution; no MPS contraction or ED; backend digest maps matched without rehashing backend files",
            "analysis_source_sha256": {label(Path(m.__file__).resolve()): digest(Path(m.__file__).resolve()) for m in modules},
            "inputs": rows, "run_verification": {k: c["verification"] for k, c in cases.items()},
            "geometry": geometry, "runtime": reference["record"]["runtime"], "stationarity_limits": LIMITS,
            "comparisons": comparisons, "field_convention": "F_Q=E_Q-h*Q/2; h_lower=E3-E1; h_upper=E5-E3; exchange energies, J=1",
            "field_intervals": intervals, "signed_sector_column_changes": columns,
            "historical_limitations": ["Q3 chi512 cumulative9 is the completed first batch of a timed-out worker; second batch is excluded",
                                        "Q3 chi256 baseline lacks raw complex second-moment components and backend_unchanged flag; absence retained"]}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--q1-diagnostics", type=Path, required=True)
    parser.add_argument("--q5-diagnostics", type=Path, required=True)
    args = parser.parse_args()
    output, temporary = args.output.resolve(), None
    try:
        require(output.suffix == ".json" and not output.exists(), "use a new .json output path")
        result = analyze(args.q1_diagnostics.resolve(), args.q5_diagnostics.resolve())
        encoded = json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
        output.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=output.parent, delete=False) as handle:
            temporary = Path(handle.name)
            handle.write(encoded)
        os.link(temporary, output)
    except (OSError, KeyError, TypeError, ValueError) as error:
        parser.exit(2, f"static chi512 follow-up: {error}\n")
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    print(json.dumps({"output": str(output), "field_intervals": result["field_intervals"]}, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
