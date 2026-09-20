#!/usr/bin/env python3
"""Compare the explicit N27 Q1/Q5 chi256 cumulative6-to8 refinement.

Only completed, verified batches enter this finite-trial comparison. No MPS
contractions are performed; file hashes and saved measurements are checked.
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
import summarize_static_campaign as campaign

ROOT = Path(__file__).resolve().parent.parent
LIMITS = {"energy_per_site": 1e-6, "sz_profile": 1e-4, "bond_profile": 1e-4,
          "variance_per_site": 1e-5, "truncation": 1e-6}
ROLES = {"old_q1": (1, 6), "old_q5": (5, 6), "fixed_q3": (3, 8),
         "new_q1": (1, 8), "new_q5": (5, 8)}


def require(condition: bool, message: str) -> None:
    if not condition:
        raise ValueError(message)


def rootpath(value: str) -> Path:
    path = Path(value)
    return (path if path.is_absolute() else ROOT / path).resolve()


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def label(path: Path) -> str:
    return os.path.relpath(path, ROOT)


def load_input(path: Path, role: str) -> dict:
    path = path.resolve()
    payload = path.read_bytes()
    sha = hashlib.sha256(payload).hexdigest()
    record = tomllib.loads(payload.decode("utf-8"))
    q, sweeps = ROLES[role]
    require(record.get("schema_version") == 1, f"{role}: unsupported schema")
    require(record.get("N") == 27 and record.get("Q") == q, f"{role}: wrong N/Q")
    geometry = campaign._geometry(record["configuration"], record)
    require((geometry["Lx"], geometry["Ly"]) == (3, 3), f"{role}: wrong geometry")
    require(geometry.get("gauge") == "seam" and
            geometry.get("ordering") == "x_then_y_then_A_B_C", f"{role}: unsupported convention")
    cfg, solver = record["config"], record["solver"]
    require(all(cfg[k] == v for k, v in {"Lx": 3, "Ly": 3, "Q": q, "maxdim": 256}.items()),
            f"{role}: configuration mismatch")
    require(cfg["stationarity"] == LIMITS, f"{role}: stationarity criteria changed")
    require(solver.get("maxdim") == 256 and solver.get("nsweeps") == 2 and
            solver.get("cutoff") == solver.get("noise") == 0, f"{role}: wrong solver")
    require(record.get("status") != "source_or_parent_changed" and
            all(record.get(k) is not False for k in ("sources_unchanged", "backend_unchanged", "parent_unchanged")),
            f"{role}: source or parent integrity failure")
    rows = record["batches"]
    eligible = [r for r in rows if campaign._eligible_batch(r)]
    require(bool(eligible), f"{role}: no fully completed verified row")
    latest = max(r["cumulative_sweeps"] for r in eligible)
    chosen = [r for r in eligible if r["cumulative_sweeps"] == latest]
    require(len(chosen) == 1 and latest == sweeps, f"{role}: required cumulative{sweeps} final row unavailable")
    row = chosen[0]
    require(row.get("maxlinkdim") == 256, f"{role}: actual chi is not 256")
    energy = campaign._number(row["energy"], "energy")
    variance = campaign._number(row["variance"], "variance")
    sz = campaign._vector(row["sz_profile"], 27, "sz_profile")
    bonds = campaign._vector(row["bond_energy"], len(geometry["bonds"]), "bond_energy")
    trunc = [campaign._number(x, "truncation") for x in row["measured_truncation_errors"]]
    require(len(trunc) == 2 and all(0 <= x <= 1 for x in trunc), f"{role}: missing/invalid truncation")
    require(abs(math.fsum(sz) - q / 2) <= 1e-10 and abs(math.fsum(bonds) - energy) <= 1e-8,
            f"{role}: profile sum rule failed")
    cp = rootpath(row["checkpoint"])
    hashes = {name: digest(cp / name) for name in ("metadata.toml", "state.jls", "checksums.toml")}
    require(hashes["metadata.toml"] == row["checkpoint_metadata_sha256"] and
            hashes["state.jls"] == row["checkpoint_payload_sha256"], f"{role}: checkpoint digest mismatch")
    backend = record["code"]["source_sha256"]
    require(isinstance(backend, dict) and bool(backend), f"{role}: backend hashes absent")
    require(all(isinstance(v, str) and len(v) == 64 and all(c in "0123456789abcdef" for c in v)
                for v in backend.values()), f"{role}: invalid backend digest")
    execution, audit = campaign._companion_provenance(path, sha)
    summary = {
        "validation_path": label(path), "validation_sha256": sha, "case_id": record["case_id"],
        "worker_status": record.get("status"), "execution": execution, "post_run_audit": audit,
        "source_finalization": {k: record.get(k) for k in ("sources_unchanged", "backend_unchanged", "parent_unchanged")},
        "batch_statuses": [{k: r.get(k) for k in ("batch", "cumulative_sweeps", "status", "active_phase",
                                                "integrity_passed", "checkpoint_verified")} for r in rows],
        "selected_batch": row["batch"], "cumulative_sweeps": sweeps,
        "checkpoint": label(cp), "checkpoint_sha256": hashes, "Q": q, "N": 27,
        "maxdim": 256, "energy": energy, "variance": variance, "variance_per_site": variance / 27,
        "last_sweep_truncation_error": trunc[-1], "sz_profile": sz, "bond_energy": bonds,
        "quality": campaign._quality_summary(row),
        "backend_source_sha256": backend,
        "recorded_analysis_source_sha256": {k: v for k, v in record.get("analysis_source_sha256", {}).items()
                                            if k in ("examples/research_static.jl", "examples/run_research.py",
                                                     "examples/run_static18.py", "examples/run_static27.py")},
    }
    return {"path": path, "record": record, "geometry": geometry, "row": row, "summary": summary}


def check_parent(old: dict, new: dict) -> dict:
    record, before = new["record"], old["summary"]
    cfg, parent = record["config"], record["parent"]
    cp = rootpath(before["checkpoint"])
    require(cfg["initialization"] == "resume" and
            cfg["initial_completed_sweeps"] == parent["initial_completed_sweeps"] == 6,
            "new case did not resume cumulative6")
    require(parent.get("strict_source_runtime_configuration_load") is True, "strict parent load absent")
    require(rootpath(cfg["parent_checkpoint"]) == rootpath(parent["checkpoint"]) == cp,
            "parent checkpoint link mismatch")
    evidence = rootpath(parent["record"])
    require(rootpath(cfg["parent_record"]) == evidence and
            digest(evidence) == parent["record_sha256"] == cfg["parent_record_sha256"] == before["validation_sha256"],
            "parent evidence digest mismatch")
    require(cfg["parent_metadata_sha256"] == before["checkpoint_sha256"]["metadata.toml"],
            "parent metadata digest mismatch")
    recorded = {rootpath(k): v for k, v in parent["checkpoint_sha256"].items()}
    expected = {cp / name: sha for name, sha in before["checkpoint_sha256"].items()}
    expected[evidence] = before["validation_sha256"]
    require(recorded == expected, "parent checkpoint/evidence digest table mismatch")
    require(abs(parent["loaded_energy"] - before["energy"]) <= 1e-12, "loaded parent energy mismatch")
    return {"status": "verified", "initial_completed_sweeps": 6,
            "checkpoint": label(cp), "evidence": label(evidence), "evidence_sha256": before["validation_sha256"]}


def compare(old: dict, new: dict) -> dict:
    a, b = old["summary"], new["summary"]
    delta_e = b["energy"] - a["energy"]
    sz = campaign._difference(a["sz_profile"], b["sz_profile"])
    bonds = campaign._difference(a["bond_energy"], b["bond_energy"])
    last_sweep = new["row"].get("last_sweep_energy_change")
    values = {"energy_per_site": None if last_sweep is None else max(abs(delta_e), abs(last_sweep)) / 27,
              "sz_profile": max(abs(x - y) for x, y in zip(a["sz_profile"], b["sz_profile"])),
              "bond_profile": max(abs(x - y) for x, y in zip(a["bond_energy"], b["bond_energy"])),
              "variance_per_site": abs(b["variance_per_site"]), "truncation": b["last_sweep_truncation_error"]}
    criteria = {k: {"value": v, "limit": LIMITS[k], "passed": None if v is None else v <= LIMITS[k]}
                for k, v in values.items()}
    all_passed = None if any(v["passed"] is None for v in criteria.values()) else all(v["passed"] for v in criteria.values())
    reported = new["row"].get("stationarity_passed")
    if all_passed is not None and reported is not None:
        require(all_passed == reported, "stationarity recomputation disagrees with worker")
    return {"Q": b["Q"], "parent_link": check_parent(old, new),
            "energy_change": delta_e, "energy_change_per_site": delta_e / 27,
            "variance_change": b["variance"] - a["variance"],
            "last_truncation_change": b["last_sweep_truncation_error"] - a["last_sweep_truncation_error"],
            "scalar_change_convention": "new minus old",
            "profile_difference_convention": "old minus new at identical site/bond coordinates; mean_difference retains this sign",
            "sz_profile_difference": sz, "bond_profile_difference": bonds,
            "stationarity_criteria": criteria, "stationarity_passed_recomputed": all_passed,
            "stationarity_passed_recorded": reported,
            "y_shift_only_comparison": campaign._circumference_shift_comparison(old["geometry"], a, b, 0.0)}


def interval(q1: dict, q3: dict, q5: dict) -> dict:
    lower, upper = q3["energy"] - q1["energy"], q5["energy"] - q3["energy"]
    return {"h_lower": lower, "h_upper": upper, "width": upper - lower,
            "finite_width_in_compared_trial_set": lower < upper}


def analyze(paths: dict[str, Path]) -> dict:
    cases = {role: load_input(paths[role], role) for role in ROLES}
    reference = cases["fixed_q3"]
    geometry = {k: v for k, v in reference["geometry"].items() if k != "Q"}
    for role, case in cases.items():
        require({k: v for k, v in case["geometry"].items() if k != "Q"} == geometry, f"{role}: geometry mismatch")
        require(case["record"]["solver"] == reference["record"]["solver"], f"{role}: solver mismatch")
        require(case["summary"]["backend_source_sha256"] == reference["summary"]["backend_source_sha256"],
                f"{role}: backend source mismatch")
        require(case["record"]["runtime"] == reference["record"]["runtime"], f"{role}: runtime mismatch")
    rows = {k: v["summary"] for k, v in cases.items()}
    before = interval(rows["old_q1"], rows["fixed_q3"], rows["old_q5"])
    after = interval(rows["new_q1"], rows["fixed_q3"], rows["new_q5"])
    columns = {}
    for stage, prefix in (("before", "old"), ("after", "new")):
        pairs = [(f"{prefix}_q1", "fixed_q3"), ("fixed_q3", f"{prefix}_q5")]
        columns[stage] = []
        for lo, hi in pairs:
            result = campaign._signed_column_delta(reference["geometry"], rows[lo]["sz_profile"], rows[hi]["sz_profile"], 2)
            require(result["sum_rule_passed"], "column charge sum rule failed")
            columns[stage].append({"lower_Q": rows[lo]["Q"], "higher_Q": rows[hi]["Q"], **result})
    return {"schema_version": 1, "recorded_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
            "study": "N27 Q1/Q5 cumulative6_to8 chi256; Q3 cumulative8 fixed",
            "claims": {"evidence": "finite variational trials", "rigorous_boundary_error_bounds": False,
                       "matched_accuracy_inferred": False, "ground_state_convergence": "not_established",
                       "bulk_plateau": "not_established", "phase_identification": "not_attempted",
                       "unexplored_sectors": "Q>=7; no exclusion of sector skipping",
                       "variance": "Hamiltonian dispersion, not a ground-energy error bound"},
            "verification_scope": "saved scalar/profile arithmetic and file hashes; no new MPS contraction or ED; backend digest maps compared but backend sources not rehashed here",
            "analysis_source_sha256": {label(Path(__file__).resolve()): digest(Path(__file__).resolve()),
                                       label(Path(campaign.__file__).resolve()): digest(Path(campaign.__file__).resolve())},
            "inputs": rows, "solver": reference["record"]["solver"], "geometry": geometry,
            "stationarity_limits": LIMITS, "field_convention": "F_Q=E_Q-h*Q/2; total exchange energies; J=1",
            "field_intervals": {"before": before, "after": after,
                                "change": {k: after[k] - before[k] for k in ("h_lower", "h_upper", "width")}},
            "refinements": [compare(cases[f"old_q{q}"], cases[f"new_q{q}"]) for q in (1, 5)],
            "signed_sector_column_changes": columns}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    for role in ROLES:
        parser.add_argument("--" + role.replace("_", "-"), type=Path, required=True)
    args = parser.parse_args()
    paths = {role: getattr(args, role) for role in ROLES}
    output = args.output.resolve()
    try:
        require(output.suffix == ".json", "output must be JSON")
        require(not output.exists(), "use a new output path")
        result = analyze(paths)
        encoded = json.dumps(result, indent=2, sort_keys=True, allow_nan=False) + "\n"
        output.parent.mkdir(parents=True, exist_ok=True)
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=output.parent, delete=False) as handle:
            temporary = Path(handle.name)
            handle.write(encoded)
        os.replace(temporary, output)
    except (OSError, KeyError, TypeError, ValueError) as error:
        parser.exit(2, f"sector refinement analysis: {error}\n")
    print(json.dumps({"output": str(output), "field_intervals": result["field_intervals"]}, allow_nan=False))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
