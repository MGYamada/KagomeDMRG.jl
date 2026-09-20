#!/usr/bin/env python3
"""Summarize atomic research_static.jl validation snapshots, using stdlib only.

Usage: python3 examples/summarize_static_campaign.py OUTPUT.toml VALIDATION.toml ...
       python3 examples/summarize_static_campaign.py --selfcheck

This is a descriptive analysis of finite-chi trials. Integrity, numerical
stationarity, and ground-state convergence are separate. Every eligible batch
is retained; the latest eligible batch represents each supplied case. No
energy minimum, field interval, or spatial profile establishes a bulk phase.
Sibling execution and post-run audit records retain their separate reported
statuses. Reading an audit does not rerun its historical source-file checks.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import itertools
import json
import math
import os
from pathlib import Path
import tempfile
import time
import tomllib


ROOT = Path(__file__).resolve().parent.parent
COMPLETE_BATCH = "diagnostics_passed_accuracy_separate"
CONFIG_KEYS = {
    "schema_version", "case_id", "Lx", "Ly", "Q", "seed", "maxdim", "batches",
    "batch_sweeps", "initialization", "parent_checkpoint", "parent_record",
    "parent_record_sha256", "parent_metadata_sha256", "initial_completed_sweeps",
    "backend_snapshot_origin", "stationarity", "helper",
}
SOLVER_KEYS = {
    "seed", "nsweeps", "maxdim", "cutoff", "noise", "eigsolve_tol",
    "eigsolve_krylovdim", "eigsolve_maxiter", "measure_variance",
}
GEOMETRY_KEYS = {
    "Lx", "Ly", "N", "Q", "charge_convention", "ordering", "axis_boundary",
    "circumference_boundary", "termination", "wrap", "exchange_phase",
    "bond_family_convention", "gauge", "hz", "sites", "bonds",
}
EXECUTION_STRING_KEYS = {
    "status", "worker_validation_status", "resource_measurement_method",
    "resource_measurement_scope", "resource_cpu_method", "resource_peak_rss_method",
    "resource_ru_maxrss_native_unit", "resource_measurement_status",
}
EXECUTION_BOOL_KEYS = {
    "worker_exit_confirmed", "termination_grace_included_in_wall_limit",
    "resource_excludes_launcher", "resource_is_simultaneous_process_tree_peak",
    "resource_prior_child_usage_detected",
}
EXECUTION_INT_KEYS = {
    "schema_version", "julia_threads", "blas_threads", "worker_exit_code",
    "worker_peak_rss_bytes", "resource_ru_maxrss_bytes_per_native_unit",
}
EXECUTION_NUMBER_KEYS = {
    "wall_limit_seconds", "termination_grace_seconds", "sigterm_elapsed_seconds",
    "sigkill_elapsed_seconds", "wall_elapsed_seconds", "worker_cpu_user_seconds",
    "worker_cpu_system_seconds", "worker_cpu_total_seconds",
}


def _digest(value: object) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                                     allow_nan=False).encode()).hexdigest()


def _integer(value: object, name: str) -> int:
    if isinstance(value, bool) or not isinstance(value, int):
        raise ValueError(f"{name} must be an integer")
    return value


def _number(value: object, name: str) -> float:
    if isinstance(value, bool) or not isinstance(value, (int, float)) or not math.isfinite(value):
        raise ValueError(f"{name} must be finite and real")
    return float(value)


def _vector(value: object, length: int, name: str) -> list[float]:
    if not isinstance(value, list) or len(value) != length:
        raise ValueError(f"{name} must have length {length}")
    return [_number(v, name) for v in value]


def _companion_paths(validation: Path) -> tuple[Path, Path]:
    """Map original output records and explicitly named archived copies only."""
    suffix = "_validation.toml"
    prefix = validation.name[:-len(suffix)] + "_" if validation.name.endswith(suffix) else ""
    return (validation.with_name(prefix + "execution.toml"),
            validation.with_name(prefix + "post_run_hash_audit.json"))


def _optional_snapshot(path: Path) -> tuple[dict, bytes | None]:
    metadata = {"path": str(path.resolve())}
    try:
        payload = path.read_bytes()
    except FileNotFoundError:
        return {**metadata, "load_status": "not_recorded"}, None
    except OSError as error:
        return {**metadata, "load_status": "unreadable", "error_type": type(error).__name__}, None
    return {**metadata, "load_status": "recorded", "sha256": hashlib.sha256(payload).hexdigest(),
            "bytes": len(payload)}, payload


def _execution_summary(record: dict) -> dict:
    """Allowlist resource/termination metadata; do not copy command or environment."""
    result = {}
    for key in EXECUTION_STRING_KEYS | EXECUTION_BOOL_KEYS | EXECUTION_INT_KEYS | EXECUTION_NUMBER_KEYS:
        if key not in record:
            continue
        value = record[key]
        if key in EXECUTION_STRING_KEYS and not isinstance(value, str):
            raise ValueError(f"execution {key} must be a string")
        if key in EXECUTION_BOOL_KEYS and not isinstance(value, bool):
            raise ValueError(f"execution {key} must be boolean")
        if key in EXECUTION_INT_KEYS:
            value = _integer(value, f"execution {key}")
        if key in EXECUTION_NUMBER_KEYS:
            value = _number(value, f"execution {key}")
        result[key] = value
    return result


def _audit_digest_link(recorded: object, actual: str | None) -> dict:
    result = {}
    if isinstance(recorded, str):
        result["recorded_sha256"] = recorded
    if actual is not None:
        result["actual_sha256"] = actual
    if recorded is None:
        status = "not_recorded"
    elif (not isinstance(recorded, str) or len(recorded) != 64
          or any(c not in "0123456789abcdefABCDEF" for c in recorded)):
        status = "invalid_recorded_digest"
    elif actual is None:
        status = "actual_file_unavailable"
    else:
        status = "matched" if recorded.lower() == actual else "mismatch"
    return {**result, "status": status}


def _post_run_audit_summary(record: dict, validation_sha256: str, execution_sha256: str | None) -> dict:
    if not isinstance(record, dict):
        raise ValueError("post-run audit must be a JSON object")
    reported_status = record.get("status", "not_recorded")
    if not isinstance(reported_status, str):
        raise ValueError("post-run audit status must be a string")
    validation_link = _audit_digest_link(record.get("validation_sha256"), validation_sha256)
    execution_link = _audit_digest_link(record.get("execution_sha256"), execution_sha256)
    checks = record.get("checks")
    valid_checks = isinstance(checks, dict) and bool(checks)
    check_values = list(checks.values()) if isinstance(checks, dict) else []
    result = {
        "reported_status": reported_status,
        "validation_digest_link": validation_link, "execution_digest_link": execution_link,
        "digest_links_match": validation_link["status"] == execution_link["status"] == "matched",
        "recorded_checks_status": "recorded" if valid_checks else "missing_empty_or_invalid",
        "recorded_check_count": len(check_values),
        "recorded_true_check_count": sum(value is True for value in check_values),
        "recorded_false_check_count": sum(value is False for value in check_values),
        "recorded_nonboolean_check_count": sum(not isinstance(value, bool) for value in check_values),
        "recorded_checks_all_true": valid_checks and all(value is True for value in check_values),
        "source_files_rehashed_by_summary": False,
        "verification_scope": "audit file hashed; recorded validation/execution digests compared with the exact input snapshots; source-file checks are historical audit reports, not rerun here",
    }
    if "worker_finalization" in record:
        if not isinstance(record["worker_finalization"], str):
            raise ValueError("audit worker_finalization must be a string")
        result["reported_worker_finalization"] = record["worker_finalization"]
    if "completed_diagnostic_sweeps" in record:
        sweeps = record["completed_diagnostic_sweeps"]
        if not isinstance(sweeps, list):
            raise ValueError("audit completed_diagnostic_sweeps must be a list")
        result["reported_completed_diagnostic_sweeps"] = [_integer(v, "audit completed sweep") for v in sweeps]
    if "source_state_checked_before_later_analysis_driver_revision" in record:
        value = record["source_state_checked_before_later_analysis_driver_revision"]
        if not isinstance(value, bool):
            raise ValueError("audit source-state flag must be boolean")
        result["reported_source_state_checked_before_later_analysis_driver_revision"] = value
    return result


def _companion_provenance(validation: Path, validation_sha256: str) -> tuple[dict, dict]:
    execution_path, audit_path = _companion_paths(validation)
    execution, payload = _optional_snapshot(execution_path)
    if payload is not None:
        try:
            execution["reported"] = _execution_summary(tomllib.loads(payload.decode("utf-8")))
        except (UnicodeError, TypeError, ValueError) as error:
            execution.update(load_status="invalid_recorded", error_type=type(error).__name__)
    audit, payload = _optional_snapshot(audit_path)
    if payload is not None:
        try:
            audit.update(_post_run_audit_summary(json.loads(payload.decode("utf-8")),
                                                validation_sha256, execution.get("sha256")))
        except (UnicodeError, TypeError, ValueError) as error:
            audit.update(load_status="invalid_recorded", error_type=type(error).__name__)
    return execution, audit


def _stats(values: list[float], weights: list[float] | None = None) -> dict:
    """RMS is uncentered; rms_about_mean is a population standard deviation."""
    if not values:
        return {"count": 0, "status": "empty_selection"}
    weights = [1.0] * len(values) if weights is None else weights
    total_weight = math.fsum(weights)
    mean = math.fsum(w * x for w, x in zip(weights, values)) / total_weight
    return {
        "count": len(values), "weight_sum": total_weight,
        "weighted_sum": math.fsum(w * x for w, x in zip(weights, values)),
        "mean": mean,
        "rms": math.sqrt(math.fsum(w * x * x for w, x in zip(weights, values)) / total_weight),
        "rms_about_mean": math.sqrt(math.fsum(w * (x - mean) ** 2 for w, x in zip(weights, values)) / total_weight),
        "minimum": min(values), "maximum": max(values),
    }


def _difference(a: list[float], b: list[float]) -> dict:
    if len(a) != len(b) or not a:
        raise ValueError("profile differences require equally sized, nonempty matched profiles")
    delta = [x - y for x, y in zip(a, b)]
    return {"rms_difference": math.sqrt(math.fsum(x * x for x in delta) / len(delta)),
            "max_abs_difference": max(map(abs, delta)), "mean_difference": math.fsum(delta) / len(delta)}


def _eligible_batch(row: dict) -> bool:
    return (row.get("status") == COMPLETE_BATCH and row.get("active_phase") == "completed"
            and row.get("integrity_passed") is True and row.get("checkpoint_verified") is True)


def _quality_summary(row: dict) -> dict:
    """Preserve missing measurements and distinguish fixed-chi comparisons.

    A one-sweep batch can have an across-batch stationarity comparison while
    its within-batch sweep change is unmeasured. Legacy records do not prove
    that a reported stationarity result compared the same configured chi.
    """
    for key in ("within_batch_sweep_change_measured", "same_chi_comparison",
                "stationarity_evaluated", "stationarity_passed"):
        if key in row and not isinstance(row[key], bool):
            raise ValueError(f"{key} must be boolean when recorded")
    has_change = "last_sweep_energy_change" in row
    measured = row.get("within_batch_sweep_change_measured", has_change)
    if measured != has_change:
        raise ValueError("within-batch sweep-change flag contradicts the recorded value")
    result = {
        "within_batch_sweep_change_measured": measured,
        "within_batch_sweep_change_flag_source": "recorded" if "within_batch_sweep_change_measured" in row else "legacy_value_presence",
        "HdaggerH_real_status": "recorded" if "HdaggerH_real" in row else "not_recorded",
        "HdaggerH_imag_status": "recorded" if "HdaggerH_imag" in row else "not_recorded",
    }
    for key in ("last_sweep_energy_change", "last_optimizer_energy_error", "HdaggerH_real", "HdaggerH_imag"):
        if key in row:
            result[key] = _number(row[key], key)
    for key in ("same_chi_comparison", "stationarity_evaluated", "stationarity_passed"):
        if key in row:
            result[f"{key}_recorded"] = row[key]
    kind = row.get("comparison_kind", "not_recorded")
    if not isinstance(kind, str):
        raise ValueError("comparison_kind must be a string when recorded")
    result["comparison_kind"] = kind
    chi_changed = kind == "chi_change" or row.get("same_chi_comparison") is False
    same_chi = not chi_changed and (kind == "same_chi" or row.get("same_chi_comparison") is True)
    evaluated, passed = row.get("stationarity_evaluated"), row.get("stationarity_passed")
    issues = []
    if (kind == "chi_change" and row.get("same_chi_comparison") is True) or (kind == "same_chi" and row.get("same_chi_comparison") is False):
        issues.append("comparison_kind_and_same_chi_flag_disagree")
    if chi_changed:
        status = "not_evaluated_chi_change"
        if evaluated is True or passed is True:
            issues.append("fixed_chi_stationarity_claimed_for_chi_change")
    elif evaluated is False:
        status = "not_evaluated"
        if passed is True:
            issues.append("stationarity_passed_without_evaluation")
    elif passed is None:
        status = "evaluated_result_not_recorded" if evaluated is True else "not_evaluated"
    elif same_chi:
        status = "passed" if passed else "failed"
    else:
        status = "legacy_reported_passed_scope_unrecorded" if passed else "legacy_reported_failed_scope_unrecorded"
    result.update({
        "stationarity_status": status,
        "fixed_chi_stationarity_passed": status == "passed" and same_chi,
        "stationarity_metadata_consistent": not issues,
        "stationarity_metadata_issues": issues,
        "stationarity_comparison_scope": "chi_change" if chi_changed else "same_chi" if same_chi else "not_recorded",
    })
    return result


def _geometry(configuration: dict, record: dict) -> dict:
    geometry = {key: configuration[key] for key in GEOMETRY_KEYS if key in configuration}
    n, q = _integer(record["N"], "N"), _integer(record["Q"], "Q")
    lx, ly = _integer(configuration["Lx"], "Lx"), _integer(configuration["Ly"], "Ly")
    if n != 3 * lx * ly or lx < 3 or ly < 3 or abs(q) > n or (n + q) % 2:
        raise ValueError("invalid kagome cylinder or charge")
    if configuration.get("N") != n or configuration.get("Q") != q:
        raise ValueError("record and configuration disagree on N/Q")
    if configuration.get("charge_convention") != "q=2Sz":
        raise ValueError("unsupported charge convention")
    if configuration.get("axis_boundary") != "open" or configuration.get("circumference_boundary") != "periodic":
        raise ValueError("unsupported boundary conditions")
    if _number(record.get("theta"), "theta") != 0 or any(_vector(configuration["hz"], n, "hz")):
        raise ValueError("summary accepts zero-flux, zero-field exchange energies only")
    sites = configuration["sites"]
    if len(sites) != n:
        raise ValueError("site count mismatch")
    for index, site in enumerate(sites, 1):
        x, remainder = divmod(index - 1, 3 * ly)
        y, sub = divmod(remainder, 3)
        if (site.get("index"), site.get("x"), site.get("y"), site.get("sublattice")) != (index, x, y, "ABC"[sub]):
            raise ValueError("sites must use canonical x,y,A,B,C ordering")
    bonds = configuration["bonds"]
    if len(bonds) != (6 * lx - 2) * ly:
        raise ValueError("unexpected NN bond count")
    for bond in bonds:
        i, j = _integer(bond["i"], "bond i"), _integer(bond["j"], "bond j")
        if not (1 <= i <= n and 1 <= j <= n and i != j):
            raise ValueError("invalid bond endpoint")
        if bond.get("family") != "J1" or bond.get("Jxy") != 1 or bond.get("Jz") != 1:
            raise ValueError("only isotropic J1=1 NN records are summarized")
    return geometry


def _spatial_summary(geometry: dict, sz: list[float], bond: list[float]) -> tuple[list, dict]:
    sites, edges, lx = geometry["sites"], geometry["bonds"], geometry["Lx"]
    site_x = [site["x"] for site in sites]
    target = geometry["Q"] / (2 * geometry["N"])
    columns = []
    for x in range(lx):
        values = [value for value, sx in zip(sz, site_x) if sx == x]
        within = [value for value, edge in zip(bond, edges)
                  if site_x[edge["i"] - 1] == site_x[edge["j"] - 1] == x]
        incident, weights = [], []
        for value, edge in zip(bond, edges):
            weight = 0.5 * ((site_x[edge["i"] - 1] == x) + (site_x[edge["j"] - 1] == x))
            if weight:
                incident.append(value)
                weights.append(weight)
        columns.append({"x": x, "sz": _stats(values),
                        "sz_rms_about_target": math.sqrt(math.fsum((v - target) ** 2 for v in values) / len(values)),
                        "nn_within_column": _stats(within), "nn_incident_half_endpoint": _stats(incident, weights)})

    center = [lx // 2] if lx % 2 else [lx // 2 - 1, lx // 2]
    selections = {"central": center, "edge": [0, lx - 1], "interior": list(range(1, lx - 1))}
    regions = {}
    for label, selected in selections.items():
        values = [value for value, x in zip(sz, site_x) if x in selected]
        within = [value for value, edge in zip(bond, edges)
                  if site_x[edge["i"] - 1] == site_x[edge["j"] - 1] and site_x[edge["i"] - 1] in selected]
        regions[label] = {"columns": selected, "sz": _stats(values), "nn_within_column": _stats(within),
                          "sz_rms_about_target": math.sqrt(math.fsum((v - target) ** 2 for v in values) / len(values))}
    contrasts = {}
    for quantity in ("sz", "nn_within_column"):
        for metric in ("mean", "rms", "rms_about_mean"):
            contrasts[f"{quantity}_{metric}_central_minus_edge"] = regions["central"][quantity][metric] - regions["edge"][quantity][metric]
    return columns, {"regions": regions, "central_edge_contrasts": contrasts,
                     "definition": "outermost columns are edges; one/two middle columns are central; NN central/edge metrics use within-column bonds",
                     "interpretation": "descriptive finite-cylinder profiles; no edge-localization or bulk claim"}


def _sz_fourier(geometry: dict, sz: list[float]) -> dict:
    """Cell-coordinate, sublattice-resolved Fourier amplitudes, not a structure factor."""
    sites = geometry["sites"]
    sublattices = sorted({site["sublattice"] for site in sites})
    selections = {sub: [i for i, site in enumerate(sites) if site["sublattice"] == sub] for sub in sublattices}
    means = {sub: math.fsum(sz[i] for i in indices) / len(indices) for sub, indices in selections.items()}
    modes = []
    period9_power, other_power = 0.0, 0.0
    for m, n in itertools.product((-1, 0, 1), repeat=2):
        amplitudes = {}
        total_real, total_imag = 0.0, 0.0
        for sub, indices in selections.items():
            angles = [-2 * math.pi * (m * sites[i]["x"] + n * sites[i]["y"]) / 3 for i in indices]
            real = math.fsum((sz[i] - means[sub]) * math.cos(angle) for i, angle in zip(indices, angles)) / len(indices)
            imag = math.fsum((sz[i] - means[sub]) * math.sin(angle) for i, angle in zip(indices, angles)) / len(indices)
            power = real * real + imag * imag
            amplitudes[sub] = {"real": real, "imag": imag, "magnitude": math.sqrt(power), "power": power}
            total_real += len(indices) * real / len(sites)
            total_imag += len(indices) * imag / len(sites)
        mean_power = math.fsum(value["power"] for value in amplitudes.values()) / len(amplitudes)
        category = "mean_subtracted_zero_mode" if (m, n) == (0, 0) else "period9_wavevector_pair" if (m, n) in ((1, -1), (-1, 1)) else "additional_3x3_grid_component"
        if category == "period9_wavevector_pair":
            period9_power += mean_power
        elif category == "additional_3x3_grid_component":
            other_power += mean_power
        modes.append({"m": m, "n": n, "k1_over_2pi": m / 3, "k2_over_2pi": n / 3, "category": category,
                      "sublattice_amplitudes": amplitudes, "sublattice_mean_power": mean_power,
                      "all_site_normalized_real": total_real, "all_site_normalized_imag": total_imag})
    return {"status": "computed", "sublattice_means_removed": means,
            "sublattice_site_counts": {sub: len(indices) for sub, indices in selections.items()},
            "normalization": "F_s(m,n)=sum_{i in s}[(Sz_i-mean_s)*exp(-2pi*i*(m*x_i+n*y_i)/3)]/N_s; no sublattice-position phase",
            "power_normalization": "mean_s(abs(F_s)^2); category powers sum this quantity over their listed modes",
            "amplitude_units": "physical Sz", "power_units": "physical Sz squared",
            "modes": modes, "period9_pair_power": period9_power, "additional_grid_power": other_power,
            "interpretation": "OBC x Fourier amplitudes describe finite profiles, not conserved momentum or bulk order; one additional mode alone does not prove a 27-site primitive cell"}


def _connected_zz_windows(geometry: dict, sz: list[float], row: dict) -> dict:
    """Uniformly weighted connected correlations in explicitly listed coordinate windows."""
    source = row.get("correlations", {})
    common = {"source": "correlations.zz_real", "definition": "Czz_ij=Re<Sz_i Sz_j>-<Sz_i><Sz_j>",
              "normalization": "uniform mean and sqrt(mean(Czz^2)) over the selected pairs; i=j excluded",
              "coordinate_rule": "x-column windows include every y and sublattice",
              "interpretation": "finite-state correlations, not a bulk correlation-length or edge-localization claim; fixed Q constrains full row sums"}
    if not isinstance(source, dict) or "zz_real" not in source:
        return {**common, "status": "not_recorded", "missing_value_policy": "no zero filling or reconstruction"}
    count = len(sz)
    try:
        matrix = source["zz_real"]
        if not isinstance(matrix, list) or len(matrix) != count:
            raise ValueError("zz_real must have N rows")
        zz = [_vector(values, count, "zz_real row") for values in matrix]
    except (TypeError, ValueError) as error:
        return {**common, "status": "invalid_recorded_matrix", "reason": str(error)}
    connected = [[zz[i][j] - sz[i] * sz[j] for j in range(count)] for i in range(count)]
    lx = geometry["Lx"]
    center = [lx // 2] if lx % 2 else [lx // 2 - 1, lx // 2]
    columns = {"left_edge": [0], "right_edge": [lx - 1], "central": center,
               "edges": [0, lx - 1], "interior": list(range(1, lx - 1))}
    windows = []
    for left, right in (("left_edge", "central"), ("central", "right_edge"), ("left_edge", "right_edge"),
                        ("central", "edges"), ("central", "central"), ("interior", "interior")):
        left_sites = [i for i, site in enumerate(geometry["sites"]) if site["x"] in columns[left]]
        right_sites = [i for i, site in enumerate(geometry["sites"]) if site["x"] in columns[right]]
        values = [connected[i][j] for i in left_sites for j in right_sites if i != j and (left != right or i < j)]
        windows.append({"left_window": left, "right_window": right, "left_columns": columns[left], "right_columns": columns[right],
                        "left_site_count": len(left_sites), "right_site_count": len(right_sites),
                        "pair_rule": "unique distinct pairs i<j" if left == right else "all cross-window pairs",
                        "statistics": _stats(values)})
    return {**common, "status": "computed", "windows": windows,
            "fixed_Q_connected_row_sum_max_abs": max(abs(math.fsum(values)) for values in connected)}


def _row_summary(row: dict, case: dict, input_id: str) -> dict:
    geometry, config = case["configuration"], case["config"]
    n, q = case["N"], case["Q"]
    sz = _vector(row["sz_profile"], n, "sz_profile")
    bonds = _vector(row["bond_energy"], len(geometry["bonds"]), "bond_energy")
    energy, variance = _number(row["energy"], "energy"), _number(row["variance"], "variance")
    truncations = row["measured_truncation_errors"]
    if not isinstance(truncations, list) or not truncations:
        raise ValueError("missing measured truncation errors")
    truncations = [_number(value, "truncation") for value in truncations]
    if not all(0 <= value <= 1 for value in truncations):
        raise ValueError("truncation error outside [0,1]")
    if abs(math.fsum(sz) - q / 2) > 1e-10 or abs(math.fsum(bonds) - energy) > 1e-8:
        raise ValueError("saved profiles fail charge or bond-energy sum check")
    sweeps = _integer(row["cumulative_sweeps"], "cumulative_sweeps")
    actual_chi = _integer(row["maxlinkdim"], "maxlinkdim")
    if sweeps <= 0 or not (0 < actual_chi <= config["maxdim"]):
        raise ValueError("invalid sweep count or actual bond dimension")
    columns, spatial = _spatial_summary(geometry, sz, bonds)
    summary = {
        "candidate_id": f"{case['case_id']}:{input_id[:12]}:sweeps{sweeps}",
        "batch": _integer(row["batch"], "batch"), "cumulative_sweeps": sweeps,
        "energy": energy, "energy_per_site": energy / n, "variance": variance,
        "variance_per_site": variance / n, "actual_maxlinkdim": actual_chi,
        "measured_truncation_errors": truncations, "max_truncation_error": max(truncations),
        "last_sweep_truncation_error": truncations[-1],
        "integrity_passed": True, "checkpoint_verified": True, "ground_state_convergence": "not_established_by_this_summary",
        "checkpoint": str(row["checkpoint"]), "sz_profile": sz, "bond_energy": bonds,
        "columns": columns, "spatial": spatial,
        "sz_fourier_3x3": _sz_fourier(geometry, sz),
        "connected_zz_windows": _connected_zz_windows(geometry, sz, row),
    }
    summary.update(_quality_summary(row))
    for key in ("energy_change_from_previous", "max_sz_change", "max_bond_change", "variance_roundoff_scale",
                "dmrg_seconds", "diagnostic_seconds", "variance_seconds"):
        if key in row:
            summary[key] = _number(row[key], key)
    for key in ("checkpoint_metadata_sha256", "checkpoint_payload_sha256"):
        if key in row:
            summary[key] = str(row[key])
    return summary


def _case_summary(record: dict, source: dict) -> dict:
    if record.get("schema_version") != 1 or not isinstance(record.get("batches"), list):
        raise ValueError("expected research_static.jl schema 1 with batches")
    if record.get("sources_unchanged") is False or record.get("parent_unchanged") is False or record.get("status") == "source_or_parent_changed":
        raise ValueError("worker recorded a source/parent integrity failure")
    config = {key: record["config"][key] for key in CONFIG_KEYS if key in record["config"]}
    geometry = _geometry(record["configuration"], record)
    for key in ("Lx", "Ly", "Q"):
        if config[key] != geometry[key]:
            raise ValueError(f"config and geometry disagree on {key}")
    for key in ("maxdim", "seed"):
        _integer(config[key], key)
    solver = {key: record["solver"][key] for key in SOLVER_KEYS if key in record["solver"]}
    if solver.get("maxdim") != config["maxdim"] or solver.get("seed") != config["seed"]:
        raise ValueError("config and solver disagree on maxdim/seed")
    case = {
        "case_id": str(record["case_id"]), "input_sha256": source["sha256"], "input_path": source["path"],
        "input_status": str(record.get("status", "unknown")), "theta": _number(record["theta"], "theta"),
        "N": geometry["N"], "Q": geometry["Q"],
        "Lx": geometry["Lx"], "Ly": geometry["Ly"], "maxdim": config["maxdim"], "seed": config["seed"],
        "initialization": str(config["initialization"]),
        "initialization_origin": str(record.get("seed_metadata", {}).get("kind", config["initialization"])),
        "config": config, "configuration": geometry, "solver": solver,
        "config_sha256_recorded": str(record.get("config_sha256", "unavailable")),
        "source_finalization": "confirmed_unchanged" if record.get("sources_unchanged") is True and record.get("parent_unchanged") is True else "not_yet_finalized",
        "backend_source_sha256": record.get("code", {}).get("source_sha256", {}),
        "analysis_source_sha256_recorded": record.get("analysis_source_sha256", {}),
        "batches": [], "excluded_batches": [],
    }
    case["physical_configuration_sha256"] = _digest(geometry)
    previous = None
    for row in sorted(record["batches"], key=lambda item: item.get("cumulative_sweeps", -1)):
        if not _eligible_batch(row):
            case["excluded_batches"].append({"batch": row.get("batch", -1), "status": str(row.get("status", "missing")),
                                             "active_phase": str(row.get("active_phase", "missing")),
                                             "reason": "not_completed_and_integrity_verified"})
            continue
        try:
            summary = _row_summary(row, case, source["sha256"])
        except (KeyError, TypeError, ValueError) as error:
            case["excluded_batches"].append({"batch": row.get("batch", -1), "reason": f"invalid_completed_row: {error}"})
            continue
        if previous is not None:
            summary["change_from_previous_selected_batch"] = {
                "reference_candidate_id": previous["candidate_id"], "energy_difference": summary["energy"] - previous["energy"],
                "sz": _difference(summary["sz_profile"], previous["sz_profile"]),
                "bond": _difference(summary["bond_energy"], previous["bond_energy"]),
            }
        case["batches"].append(summary)
        previous = summary
    case["selection_status"] = "latest_completed_integrity_verified_batch" if case["batches"] else "no_eligible_completed_batch"
    if case["batches"]:
        case["selected_candidate_id"] = case["batches"][-1]["candidate_id"]
    if case["initialization_origin"] == "resume":
        case["initialization_origin_note"] = "original initialization is not inferred from a resumed seed integer"
    return case


def _candidate(case: dict) -> dict:
    row = case["batches"][-1]
    result = {"candidate_id": row["candidate_id"], "case_id": case["case_id"], "N": case["N"], "Q": case["Q"],
            "maxdim": case["maxdim"], "actual_maxlinkdim": row["actual_maxlinkdim"], "seed": case["seed"],
            "initialization": case["initialization"], "initialization_origin": case["initialization_origin"],
            "cumulative_sweeps": row["cumulative_sweeps"], "energy": row["energy"],
            "energy_per_site": row["energy_per_site"], "variance_per_site": row["variance_per_site"],
            "max_truncation_error": row["max_truncation_error"], "stationarity_status": row["stationarity_status"],
            "source_finalization": case["source_finalization"]}
    if "input_status" in case:
        result["input_status"] = case["input_status"]
    if "execution" in case:
        result["execution_load_status"] = case["execution"]["load_status"]
        reported = case["execution"].get("reported", {})
        for key in ("status", "worker_exit_confirmed"):
            if key in reported:
                result[f"execution_reported_{key}"] = reported[key]
    for key in ("fixed_chi_stationarity_passed", "stationarity_comparison_scope", "comparison_kind",
                "same_chi_comparison_recorded", "stationarity_evaluated_recorded", "stationarity_passed_recorded",
                "stationarity_metadata_consistent", "stationarity_metadata_issues",
                "within_batch_sweep_change_measured", "HdaggerH_imag_status"):
        if key in row:
            result[key] = row[key]
    return result


def _spreads(cases: list[dict]) -> list[dict]:
    groups = {}
    for case in cases:
        groups.setdefault(case["physical_configuration_sha256"], []).append(_candidate(case))
    results = []
    for key, candidates in groups.items():
        energies = [candidate["energy"] for candidate in candidates]
        results.append({"physical_configuration_sha256": key, "N": candidates[0]["N"], "Q": candidates[0]["Q"],
                        "candidates": candidates, "candidate_count": len(candidates), "minimum_trial_energy": min(energies),
                        "maximum_trial_energy": max(energies), "energy_spread": max(energies) - min(energies),
                        "scope": "latest completed batch from each supplied case; heterogeneous quality retained",
                        "convergence": "not_inferred_from_minimum_or_spread"})
    return results


def _field_intervals(cases: list[dict]) -> list[dict]:
    groups = {}
    for case in cases:
        if (case["N"], case["Lx"], case["Ly"]) != (27, 3, 3) or case["Q"] not in (1, 3, 5):
            continue
        without_charge = {key: value for key, value in case["configuration"].items() if key != "Q"}
        groups.setdefault(_digest(without_charge), []).append(case)
    results = []
    for geometry_key, group in groups.items():
        chis = sorted({case["maxdim"] for case in group})
        selections = [(f"fixed_configured_maxdim_{chi}", [case for case in group if case["maxdim"] == chi]) for chi in chis]
        if len(chis) > 1:
            selections.append(("mixed_maxdim_descriptive_trials", group))
        for scope, selected in selections:
            sectors = {q: [_candidate(case) for case in selected if case["Q"] == q] for q in (1, 3, 5)}
            result = {"N": 27, "target_Q": 3, "geometry_model_sha256_without_Q": geometry_key,
                      "scope": scope, "configured_maxdims": sorted({case["maxdim"] for case in selected}),
                      "compared_Q": [q for q, values in sectors.items() if values],
                      "missing_Q": [q for q, values in sectors.items() if not values],
                      "sector_candidates": {str(q): values for q, values in sectors.items()},
                      "interpretation": "differences of finite-chi variational trial energies; no rigorous interval error bound or bulk plateau claim",
                      "other_magnetization_sectors": "not_tested_by_this_Q1_Q3_Q5_comparison"}
            if result["missing_Q"]:
                result["status"] = "unavailable_missing_sectors"
            else:
                lo = {q: min(c["energy"] for c in sectors[q]) for q in sectors}
                hi = {q: max(c["energy"] for c in sectors[q]) for q in sectors}
                lower, upper = lo[3] - lo[1], lo[5] - lo[3]
                result.update({"status": "computed_descriptive_trial_differences", "h_lower": lower, "h_upper": upper,
                               "interval_width": upper - lower, "selected_trial_interval_nonempty": lower < upper,
                               "h_lower_candidate_range": [lo[3] - hi[1], hi[3] - lo[1]],
                               "h_upper_candidate_range": [lo[5] - hi[3], hi[5] - lo[3]],
                               "width_candidate_range": [lo[5] + lo[1] - 2 * hi[3], hi[5] + hi[1] - 2 * lo[3]],
                               "sector_energy_spreads": {str(q): hi[q] - lo[q] for q in sectors},
                               "minimum_energy_candidate_ids": {str(q): [c["candidate_id"] for c in sectors[q] if c["energy"] == lo[q]] for q in sectors},
                               "selection": "minimum among supplied latest trials in each Q; not a convergence decision"})
            results.append(result)
    return results


def _circumference_shift_comparison(geometry: dict, left: dict, right: dict, theta: float) -> dict:
    """Compare left(x,y,s) with right(x,y+shift mod Ly,s) at fixed x.

    At zero flux the isotropic NN bond observable is symmetric in its two
    endpoints, so unordered endpoint pairs define its translation uniquely.
    Signed winding is deliberately not part of this zero-flux matching key.
    No inferred site ordering, bond sorting, or open-axis translation is used.
    """
    result = {
        "status": "unsupported_mapping", "x_shifts_applied": False,
        "shift_convention": "left(x,y,s) minus right(x,(y+shift) mod Ly,s); translate both endpoints of each bond",
        "interpretation": "circumference domain-position comparison of finite trials; not full VBC template matching, an energy/phase selection criterion, or a convergence test",
        "shifts": [],
    }
    if (theta != 0 or geometry.get("axis_boundary") != "open"
            or geometry.get("circumference_boundary") != "periodic"
            or any(geometry.get("hz", []))
            or any(b.get("family") != "J1" or b.get("Jxy") != 1 or b.get("Jz") != 1
                   for b in geometry.get("bonds", []))):
        return {**result, "status": "not_applicable_hamiltonian_or_boundary"}
    try:
        sites, bonds = geometry["sites"], geometry["bonds"]
        n, ly = len(sites), _integer(geometry["Ly"], "Ly")
        if not n or not bonds or ly < 1:
            raise ValueError("nonempty site/bond lists and positive Ly required")
        coordinates, site_ids = {}, {}
        for position, site in enumerate(sites):
            key = (_integer(site["x"], "site x"), _integer(site["y"], "site y"), site["sublattice"])
            index = _integer(site["index"], "site index")
            if key in coordinates or index in site_ids or not 0 <= key[1] < ly:
                raise ValueError("site coordinates/indices do not define a one-to-one map")
            coordinates[key], site_ids[index] = position, position
        edges = {}
        for position, bond in enumerate(bonds):
            i, j = _integer(bond["i"], "bond i"), _integer(bond["j"], "bond j")
            if i not in site_ids or j not in site_ids or i == j:
                raise ValueError("bond endpoints do not identify two sites")
            key = tuple(sorted((i, j)))
            if key in edges:
                raise ValueError("bond endpoints do not define a one-to-one map")
            edges[key] = position
        left_sz, right_sz = (_vector(row["sz_profile"], n, "sz_profile") for row in (left, right))
        left_bond, right_bond = (_vector(row["bond_energy"], len(bonds), "bond_energy") for row in (left, right))
        shifts = []
        for shift in range(ly):
            site_map = [coordinates[(site["x"], (site["y"] + shift) % ly, site["sublattice"])] for site in sites]
            if set(site_map) != set(range(n)):
                raise ValueError("translated site mapping is not a permutation")
            translated_ids = {site["index"]: sites[site_map[position]]["index"] for position, site in enumerate(sites)}
            bond_map = [edges[tuple(sorted((translated_ids[b["i"]], translated_ids[b["j"]])))] for b in bonds]
            if set(bond_map) != set(range(len(bonds))):
                raise ValueError("translated bond mapping is not a permutation")
            shifts.append({"y_shift": shift, "site_mapping_one_to_one": True, "bond_mapping_one_to_one": True,
                           "sz": _difference(left_sz, [right_sz[i] for i in site_map]),
                           "bond": _difference(left_bond, [right_bond[i] for i in bond_map])})
        minima = {}
        for quantity in ("sz", "bond"):
            minimum = min(row[quantity]["rms_difference"] for row in shifts)
            minima[quantity] = {"rms_difference": minimum,
                                "y_shifts": [row["y_shift"] for row in shifts if row[quantity]["rms_difference"] == minimum]}
        return {**result, "status": "computed", "shifts": shifts, "minimum_rms": minima,
                "site_count": n, "bond_count": len(bonds), "minimum_tie_convention": "all exact floating-point minima retained separately for Sz and bond"}
    except (KeyError, TypeError, ValueError) as error:
        return {**result, "reason": f"one-to-one coordinate/bond translation unavailable: {type(error).__name__}: {error}"}


def _profile_comparisons(cases: list[dict]) -> tuple[list, list]:
    results, excluded = [], []
    kinds = {"random", "period9", "period27"}
    for left, right in itertools.combinations(cases, 2):
        a, b = _candidate(left), _candidate(right)
        if (a["N"], a["Q"]) != (b["N"], b["Q"]) or a["initialization_origin"] == b["initialization_origin"]:
            continue
        if not {a["initialization_origin"], b["initialization_origin"]} <= kinds:
            continue
        reasons = []
        if left["physical_configuration_sha256"] != right["physical_configuration_sha256"]:
            reasons.append("Hamiltonian_or_geometry_differs")
        if left["solver"] != right["solver"]:
            reasons.append("solver_configuration_or_seed_differs")
        if reasons:
            excluded.append({"left": a["candidate_id"], "right": b["candidate_id"], "reasons": reasons})
            continue
        lrow, rrow = left["batches"][-1], right["batches"][-1]
        results.append({"left": a, "right": b, "sz": _difference(lrow["sz_profile"], rrow["sz_profile"]),
                        "bond": _difference(lrow["bond_energy"], rrow["bond_energy"]),
                        "circumference_shift_comparison": _circumference_shift_comparison(left["configuration"], lrow, rrow, left["theta"]),
                        "same_cumulative_sweeps": a["cumulative_sweeps"] == b["cumulative_sweeps"],
                        "same_backend_source_hashes": left["backend_source_sha256"] == right["backend_source_sha256"],
                        "scope": "same N,Q,Hamiltonian,geometry,solver and seed; differing completed sweep counts remain explicit",
                        "interpretation": "raw finite-trial profile difference, not an order or convergence criterion"})
    return results, excluded


def _length_comparisons(cases: list[dict]) -> list[dict]:
    results = []
    for short, long in itertools.product([case for case in cases if case["N"] == 27],
                                         [case for case in cases if case["N"] == 54]):
        if short["Q"] * long["N"] != long["Q"] * short["N"]:
            continue
        same_geometry = all(short["configuration"].get(key) == long["configuration"].get(key)
                            for key in ("Ly", "wrap", "ordering", "axis_boundary", "circumference_boundary", "termination", "gauge"))
        if not same_geometry:
            continue
        a, b = _candidate(short), _candidate(long)
        short_spatial, long_spatial = short["batches"][-1]["spatial"], long["batches"][-1]["spatial"]
        changes = {}
        for region in ("central", "edge", "interior"):
            changes[region] = {}
            for quantity in ("sz", "nn_within_column"):
                changes[region][quantity] = {f"{metric}_N54_minus_N27": long_spatial["regions"][region][quantity][metric] - short_spatial["regions"][region][quantity][metric]
                                            for metric in ("mean", "rms", "rms_about_mean")}
        results.append({"N27_candidate": a, "N54_candidate": b, "N27_spatial": short_spatial, "N54_spatial": long_spatial,
                        "aggregate_changes": changes, "same_solver_configuration": short["solver"] == long["solver"],
                        "same_initialization_origin": a["initialization_origin"] == b["initialization_origin"],
                        "same_seed": a["seed"] == b["seed"], "same_cumulative_sweeps": a["cumulative_sweeps"] == b["cumulative_sweeps"],
                        "scope": "same wrap, termination and magnetization fraction; compare coordinate-defined column aggregates only",
                        "site_index_difference_used": False,
                        "interpretation": "no bulk, VBC, or edge-localization inference; N27 has only one interior column"})
    return results


def _signed_column_delta(geometry: dict, lower: list[float], higher: list[float], delta_q: int) -> dict:
    if len(lower) != len(higher) or len(lower) != len(geometry["sites"]):
        raise ValueError("sector profile changes require the same explicit site geometry")
    columns = []
    for x in range(geometry["Lx"]):
        indices = [i for i, site in enumerate(geometry["sites"]) if site["x"] == x]
        delta = math.fsum(higher[i] - lower[i] for i in indices)
        columns.append({"x": x, "site_count": len(indices), "signed_delta_sz": delta,
                        "mean_signed_delta_sz_per_site": delta / len(indices)})
    total = math.fsum(column["signed_delta_sz"] for column in columns)
    expected = delta_q / 2
    error = abs(total - expected)
    return {"columns": columns, "signed_delta_sz_total": total, "expected_delta_sz_total": expected,
            "sum_rule_abs_error": error, "sum_rule_tolerance": 2e-10, "sum_rule_passed": error <= 2e-10,
            "normalization": "column sums of Sz(higher Q)-Sz(lower Q); no clipping, absolute value, or renormalization",
            "interpretation": "signed finite-trial charge redistribution; no edge-localization or bulk plateau inference"}


def _sector_column_changes(cases: list[dict]) -> dict:
    comparisons, availability = [], []
    for low_q, high_q in ((1, 3), (3, 5)):
        lows = [case for case in cases if case["N"] == 27 and case["Q"] == low_q]
        highs = [case for case in cases if case["N"] == 27 and case["Q"] == high_q]
        matched = 0
        for lower, higher in itertools.product(lows, highs):
            low_geometry = {key: value for key, value in lower["configuration"].items() if key != "Q"}
            high_geometry = {key: value for key, value in higher["configuration"].items() if key != "Q"}
            if low_geometry != high_geometry:
                continue
            matched += 1
            a, b = _candidate(lower), _candidate(higher)
            delta = _signed_column_delta(lower["configuration"], lower["batches"][-1]["sz_profile"],
                                         higher["batches"][-1]["sz_profile"], high_q - low_q)
            comparisons.append({"from_Q": low_q, "to_Q": high_q, "lower_Q_candidate": a, "higher_Q_candidate": b,
                                "same_solver_configuration": lower["solver"] == higher["solver"],
                                "same_initialization_origin": a["initialization_origin"] == b["initialization_origin"],
                                "same_cumulative_sweeps": a["cumulative_sweeps"] == b["cumulative_sweeps"],
                                "same_backend_source_hashes": lower["backend_source_sha256"] == higher["backend_source_sha256"],
                                "delta": delta})
        availability.append({"from_Q": low_q, "to_Q": high_q, "lower_Q_candidate_count": len(lows),
                             "higher_Q_candidate_count": len(highs), "matched_geometry_pairs": matched,
                             "status": "computed" if matched else "no_matched_eligible_sector_pair"})
    return {"status": "computed" if comparisons else "no_matched_eligible_sector_pair", "availability": availability,
            "selection": "all matched latest eligible candidate pairs, without selecting the lowest energy; chi and sweep quality remain explicit",
            "charge_convention": "Q=2M; Q1->Q3 and Q3->Q5 each require signed total DeltaSz=1",
            "comparisons": comparisons}


def summarize(paths: list[Path]) -> dict:
    inputs, cases, seen = [], [], set()
    for path in paths:
        entry = {"path": str(path.resolve())}
        try:
            # One immutable byte snapshot: hash and parse exactly the same
            # content even if the worker atomically replaces the file later.
            payload = path.read_bytes()
            entry["sha256"] = hashlib.sha256(payload).hexdigest()
            entry["bytes"] = len(payload)
            execution, audit = _companion_provenance(path, entry["sha256"])
            entry["execution"], entry["post_run_hash_audit"] = execution, audit
            if entry["sha256"] in seen:
                entry["status"] = "duplicate_input_content_not_counted_twice"
            else:
                record = tomllib.loads(payload.decode("utf-8"))
                case = _case_summary(record, entry)
                # Supplemental historical evidence never changes the original
                # validation status, source finalization, or batch eligibility.
                case["execution"], case["post_run_hash_audit"] = execution, audit
                seen.add(entry["sha256"])
                cases.append(case)
                entry["status"] = case["selection_status"]
                entry["eligible_batches"] = len(case["batches"])
                entry["excluded_batches"] = len(case["excluded_batches"])
        except (OSError, UnicodeError, tomllib.TOMLDecodeError, KeyError, TypeError, ValueError) as error:
            entry["status"] = "unavailable_or_rejected_input"
            entry["reason"] = str(error)
        inputs.append(entry)
    selected = [case for case in cases if case["batches"]]
    comparisons, excluded = _profile_comparisons(selected)
    return {
        "schema_version": 1, "recorded_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
        "status": "descriptive_summary_completed" if selected else "no_eligible_completed_batches",
        "analysis_source": "examples/summarize_static_campaign.py",
        "analysis_source_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "input_count": len(inputs), "eligible_case_count": len(selected),
        "eligible_batch_count": sum(len(case["batches"]) for case in cases),
        "charge_and_field_convention": "Q=2M, M=Q/2; F_Q(h)=E_Q-h*Q/2; h_lower=E3-E1, h_upper=E5-E3 for Q1,3,5",
        "energy_convention": "total NN isotropic J1=1 exchange energy at theta=0, hz=0",
        "candidate_policy": "retain every completed integrity-passed batch; select latest batch per supplied case; minima and spreads are descriptive only",
        "rms_convention": "rms=sqrt(mean(value^2)); rms_about_mean=sqrt(mean((value-mean)^2))",
        "measurement_coverage": "unrecorded within-batch energy changes and HdaggerH components remain absent; legacy stationarity scope is not inferred",
        "execution_provenance_policy": "execution and post-run audit are separate optional snapshots; reported timeout/running states and original source finalization are preserved; an audit report never promotes a worker to normal completion",
        "audit_verification_policy": "only the audit and its validation/execution file links are hashed here; the audit's historical source/checkpoint checks are not repeated",
        "column_bond_convention": "within-column bonds are primary; incident metrics allocate half a bond to each endpoint column",
        "claims": {"ground_state_convergence": "not_established", "bulk_plateau": "not_established",
                   "VBC": "not_identified", "edge_localization": "not_inferred", "field_error_bounds": "not_rigorous"},
        "inputs": inputs, "cases": cases, "candidate_spreads": _spreads(selected),
        "field_intervals": _field_intervals(selected), "matched_profile_comparisons": comparisons,
        "excluded_profile_comparisons": excluded, "length_comparisons": _length_comparisons(selected),
        "sector_column_magnetization_changes": _sector_column_changes(selected),
    }


def _toml_value(value: object) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, str):
        return json.dumps(value, ensure_ascii=False)
    if isinstance(value, int):
        return str(value)
    if isinstance(value, float) and math.isfinite(value):
        return repr(value)
    if isinstance(value, list):
        return "[" + ", ".join(_toml_value(item) for item in value) + "]"
    if isinstance(value, dict):
        return "{" + ", ".join(f"{_toml_value(str(key))} = {_toml_value(item)}" for key, item in sorted(value.items())) + "}"
    raise ValueError(f"unsupported or nonfinite TOML value type: {type(value).__name__}")


def dumps_toml(data: dict) -> str:
    lines = []

    def visit(table: dict, path: tuple[str, ...]) -> None:
        for key, value in table.items():
            if not isinstance(value, dict) and not (isinstance(value, list) and value and all(isinstance(item, dict) for item in value)):
                lines.append(f"{_toml_value(str(key))} = {_toml_value(value)}")
        for key, value in table.items():
            child = path + (str(key),)
            header = ".".join(_toml_value(part) for part in child)
            if isinstance(value, dict):
                lines.extend(["", f"[{header}]"])
                visit(value, child)
            elif isinstance(value, list) and value and all(isinstance(item, dict) for item in value):
                for item in value:
                    lines.extend(["", f"[[{header}]]"])
                    visit(item, child)

    visit(data, ())
    result = "\n".join(lines) + "\n"
    if tomllib.loads(result) != data:
        raise ValueError("TOML serialization round-trip failed")
    return result


def selfcheck() -> dict:
    start = time.monotonic()
    execution_path, audit_path = _companion_paths(Path("outputs/toy/validation.toml"))
    assert execution_path == Path("outputs/toy/execution.toml")
    assert audit_path == Path("outputs/toy/post_run_hash_audit.json")
    execution_path, audit_path = _companion_paths(Path("docs/research/data/toy_validation.toml"))
    assert execution_path == Path("docs/research/data/toy_execution.toml")
    assert audit_path == Path("docs/research/data/toy_post_run_hash_audit.json")
    execution = _execution_summary({"status": "timed_out", "worker_exit_confirmed": True,
                                    "wall_elapsed_seconds": 9.0, "worker_peak_rss_bytes": 1024,
                                    "command": ["omitted"], "environment": {"omitted": "omitted"}})
    assert execution == {"status": "timed_out", "worker_exit_confirmed": True,
                         "wall_elapsed_seconds": 9.0, "worker_peak_rss_bytes": 1024}
    validation_sha = hashlib.sha256(b"validation snapshot").hexdigest()
    execution_sha = hashlib.sha256(b"execution snapshot").hexdigest()
    audit_record = {"status": "passed", "validation_sha256": validation_sha,
                    "execution_sha256": execution_sha, "checks": {"source": True}}
    linked = _post_run_audit_summary(audit_record, validation_sha, execution_sha)
    assert linked["digest_links_match"] and linked["recorded_checks_all_true"]
    assert linked["reported_status"] == "passed" and not linked["source_files_rehashed_by_summary"]
    mismatch = _post_run_audit_summary(audit_record, hashlib.sha256(b"changed").hexdigest(), execution_sha)
    assert mismatch["reported_status"] == "passed" and not mismatch["digest_links_match"]
    assert mismatch["validation_digest_link"]["status"] == "mismatch"
    assert mismatch["execution_digest_link"]["status"] == "matched"
    changed_execution = _post_run_audit_summary(audit_record, validation_sha, hashlib.sha256(b"changed execution").hexdigest())
    assert changed_execution["execution_digest_link"]["status"] == "mismatch" and not changed_execution["digest_links_match"]
    missing_execution = _post_run_audit_summary(audit_record, validation_sha, None)
    assert missing_execution["execution_digest_link"]["status"] == "actual_file_unavailable"
    assert not missing_execution["digest_links_match"]
    for checks in ({}, {"source": False}, {"source": 1}):
        assert not _post_run_audit_summary({**audit_record, "checks": checks}, validation_sha, execution_sha)["recorded_checks_all_true"]
    assert tomllib.loads(dumps_toml({"execution": execution, "audit": linked, "mismatch": mismatch}))["mismatch"] == mismatch
    # Explicit three-cell circumference, with four bonds per y cell. A known
    # rightward data rotation is undone only by comparison shift +1; the bond
    # crossing y=0 verifies matching by translated endpoints across the seam.
    ring_sites = [{"index": 3 * y + sub + 1, "x": 0, "y": y, "sublattice": "ABC"[sub]}
                  for y in range(3) for sub in range(3)]
    ring_bonds = [{"i": i, "j": j, "family": "J1", "Jxy": 1, "Jz": 1}
                  for y in range(3) for i, j in ((3*y+1, 3*y+2), (3*y+1, 3*y+3),
                                                 (3*y+2, 3*y+3), (3*y+1, 3*((y-1) % 3)+3))]
    ring = {"Ly": 3, "axis_boundary": "open", "circumference_boundary": "periodic",
            "sites": ring_sites, "bonds": ring_bonds}
    ring_left = {"sz_profile": [0.1, -0.2, 0.3, -0.4, 0.05, -0.1, 0.2, 0.15, -0.25],
                 "bond_energy": [-0.3 + 0.01 * i for i in range(12)]}
    ring_right = {"sz_profile": ring_left["sz_profile"][-3:] + ring_left["sz_profile"][:-3],
                  "bond_energy": ring_left["bond_energy"][-4:] + ring_left["bond_energy"][:-4]}
    shifted = _circumference_shift_comparison(ring, ring_left, ring_right, 0.0)
    assert shifted["status"] == "computed" and not shifted["x_shifts_applied"]
    for quantity in ("sz", "bond"):
        assert shifted["minimum_rms"][quantity] == {"rms_difference": 0.0, "y_shifts": [1]}
        assert shifted["shifts"][0][quantity]["rms_difference"] > 0
        assert shifted["shifts"][1][quantity]["max_abs_difference"] == 0
    assert _circumference_shift_comparison(ring, ring_left, ring_right, 0.37)["status"] == "not_applicable_hamiltonian_or_boundary"
    assert _circumference_shift_comparison({**ring, "bonds": ring_bonds + [ring_bonds[0]]}, ring_left, ring_right, 0.0)["status"] == "unsupported_mapping"
    stats = _stats([1.0, 3.0])
    assert stats["mean"] == 2 and stats["rms_about_mean"] == 1
    assert math.isclose(stats["rms"], math.sqrt(5))
    assert _stats([1.0, 3.0], [0.5, 1.5])["mean"] == 2.5
    assert _difference([1.0, 3.0], [0.0, 1.0])["rms_difference"] == math.sqrt(2.5)
    complete = {"status": COMPLETE_BATCH, "active_phase": "completed",
                "integrity_passed": True, "checkpoint_verified": True}
    assert _eligible_batch(complete)
    for key, value in (("status", "checkpoint_saved_diagnostics_pending"), ("active_phase", "variance"),
                       ("integrity_passed", False), ("checkpoint_verified", False)):
        assert not _eligible_batch({**complete, key: value})

    def artificial_case(q: int, energy: float, chi: int = 256) -> dict:
        candidate = {"candidate_id": f"Q{q}-chi{chi}", "cumulative_sweeps": 4,
                     "energy": energy, "energy_per_site": energy / 27, "variance_per_site": 0.001,
                     "max_truncation_error": 1e-5, "actual_maxlinkdim": chi,
                     "stationarity_status": "failed"}
        return {"case_id": f"Q{q}", "N": 27, "Q": q, "Lx": 3, "Ly": 3, "maxdim": chi,
                "seed": 11, "initialization": "random", "initialization_origin": "random",
                "source_finalization": "confirmed_unchanged", "configuration": {"N": 27, "Q": q, "Lx": 3, "Ly": 3},
                "batches": [candidate]}

    field = _field_intervals([artificial_case(1, -10), artificial_case(3, -9.7), artificial_case(5, -9.2)])[0]
    assert math.isclose(field["h_lower"], 0.3, abs_tol=1e-14)
    assert math.isclose(field["h_upper"], 0.5, abs_tol=1e-14)
    assert math.isclose(field["interval_width"], 0.2, abs_tol=1e-14)
    assert _field_intervals([artificial_case(3, -9.7)])[0]["missing_Q"] == [1, 5]
    varied = _field_intervals([artificial_case(1, -10), artificial_case(3, -9.7),
                              artificial_case(3, -9.6), artificial_case(5, -9.2)])[0]
    assert all(math.isclose(x, y, abs_tol=1e-14) for x, y in zip(varied["width_candidate_range"], [0.0, 0.2]))
    assert varied["sector_candidates"]["3"][0]["stationarity_status"] == "failed"
    # A tiny explicit three-column geometry checks coordinate selections and
    # endpoint weighting without invoking any numerical solver.
    geometry = {"N": 6, "Q": 2, "Lx": 3, "sites": [{"x": x, "y": 0, "sublattice": sub} for x in range(3) for sub in "AB"],
                "bonds": [{"i": 1, "j": 2}, {"i": 3, "j": 4}, {"i": 5, "j": 6},
                          {"i": 2, "j": 3}, {"i": 4, "j": 5}]}
    columns, spatial = _spatial_summary(geometry, [0.0, 0.0, 0.5, 0.5, 0.0, 0.0], [1.0, 2.0, 3.0, 4.0, 5.0])
    assert math.isclose(sum(column["nn_incident_half_endpoint"]["weighted_sum"] for column in columns), 15.0)
    assert spatial["regions"]["central"]["columns"] == [1]
    assert spatial["regions"]["central"]["sz"]["mean"] == 0.5
    one_sweep = {**complete, "batch": 1, "cumulative_sweeps": 1, "maxlinkdim": 2,
                 "energy": 15.0, "variance": 0.01, "sz_profile": [0.0, 0.0, 0.5, 0.5, 0.0, 0.0],
                 "bond_energy": [1.0, 2.0, 3.0, 4.0, 5.0], "measured_truncation_errors": [1e-6],
                 "checkpoint": "synthetic_calibration_only", "within_batch_sweep_change_measured": False,
                 "last_optimizer_energy_error": 0.001}
    toy_case = {"case_id": "synthetic_one_sweep", "N": 6, "Q": 2,
                "configuration": geometry, "config": {"maxdim": 2}}
    one_summary = _row_summary(one_sweep, toy_case, "0" * 64)
    assert one_summary["within_batch_sweep_change_measured"] is False
    assert "last_sweep_energy_change" not in one_summary
    assert "HdaggerH_imag" not in one_summary and one_summary["HdaggerH_imag_status"] == "not_recorded"
    assert one_summary["last_optimizer_energy_error"] == 0.001
    assert one_summary["connected_zz_windows"]["status"] == "not_recorded"
    same_chi = _quality_summary({**one_sweep, "same_chi_comparison": True, "comparison_kind": "same_chi",
                                "stationarity_evaluated": True, "stationarity_passed": True})
    assert same_chi["stationarity_status"] == "passed" and same_chi["fixed_chi_stationarity_passed"]
    changed_chi = _quality_summary({**one_sweep, "same_chi_comparison": False, "comparison_kind": "chi_change",
                                   "stationarity_evaluated": False, "stationarity_passed": False})
    assert changed_chi["stationarity_status"] == "not_evaluated_chi_change"
    assert not changed_chi["fixed_chi_stationarity_passed"] and changed_chi["stationarity_metadata_consistent"]
    contradictory = _quality_summary({**one_sweep, "same_chi_comparison": False, "comparison_kind": "chi_change",
                                      "stationarity_evaluated": True, "stationarity_passed": True})
    assert contradictory["stationarity_status"] == "not_evaluated_chi_change"
    assert not contradictory["fixed_chi_stationarity_passed"] and not contradictory["stationarity_metadata_consistent"]
    unevaluated = _quality_summary({**one_sweep, "same_chi_comparison": True,
                                   "stationarity_evaluated": False, "stationarity_passed": False})
    assert unevaluated["stationarity_status"] == "not_evaluated"
    legacy = _quality_summary({"last_sweep_energy_change": 0.01, "stationarity_passed": True})
    assert legacy["stationarity_status"] == "legacy_reported_passed_scope_unrecorded"
    assert not legacy["fixed_chi_stationarity_passed"] and "HdaggerH_imag" not in legacy
    measured_moment = _quality_summary({**one_sweep, "HdaggerH_real": 225.01, "HdaggerH_imag": -1e-14})
    assert measured_moment["HdaggerH_imag"] == -1e-14 and measured_moment["HdaggerH_imag_status"] == "recorded"
    assert tomllib.loads(dumps_toml({"one_sweep": one_summary, "chi_change": changed_chi}))["one_sweep"] == one_summary
    # Values of cos/sin on the three roots of unity give hand-computable
    # Fourier coefficients. Independent (x-y) and (x+y) waves exercise the
    # period9 pair, additional grid components, normalization, and phase sign.
    fourier_geometry = {"N": 27, "Lx": 3, "Ly": 3,
                        "sites": [{"x": x, "y": y, "sublattice": sub}
                                  for x in range(3) for y in range(3) for sub in "ABC"]}
    offsets = {"A": 0.1, "B": 0.2, "C": -0.1}
    a = {"A": 0.2, "B": -0.1, "C": 0.0}
    b = {"A": 0.04, "B": 0.06, "C": -0.08}
    c = {"A": 0.06, "B": 0.0, "C": -0.02}
    cosine_values, sine_values = [1.0, -0.5, -0.5], [0.0, math.sqrt(3) / 2, -math.sqrt(3) / 2]
    pattern = [offsets[s["sublattice"]] + a[s["sublattice"]] * cosine_values[(s["x"] - s["y"]) % 3]
               + c[s["sublattice"]] * sine_values[(s["x"] - s["y"]) % 3]
               + b[s["sublattice"]] * cosine_values[(s["x"] + s["y"]) % 3] for s in fourier_geometry["sites"]]
    fourier = _sz_fourier(fourier_geometry, pattern)
    for sub in "ABC":
        assert math.isclose(fourier["sublattice_means_removed"][sub], offsets[sub], abs_tol=2e-15)
    for mode in fourier["modes"]:
        wavevector = (mode["m"], mode["n"])
        for sub, amplitude in mode["sublattice_amplitudes"].items():
            real = a[sub] / 2 if wavevector in ((1, -1), (-1, 1)) else b[sub] / 2 if wavevector in ((1, 1), (-1, -1)) else 0.0
            imag = -c[sub] / 2 if wavevector == (1, -1) else c[sub] / 2 if wavevector == (-1, 1) else 0.0
            assert math.isclose(amplitude["real"], real, abs_tol=2e-15)
            assert math.isclose(amplitude["imag"], imag, abs_tol=2e-15)
    assert math.isclose(fourier["period9_pair_power"], sum(a[sub] ** 2 + c[sub] ** 2 for sub in "ABC") / 6, abs_tol=2e-15)
    assert math.isclose(fourier["additional_grid_power"], sum(value ** 2 for value in b.values()) / 6, abs_tol=2e-15)

    # Two polarized left spins and two singlets joining center to right:
    # center/right has four pairs, two with Czz=-1/4 and two with Czz=0.
    polarized = [0.5, 0.5, 0.0, 0.0, 0.0, 0.0]
    zz = [[0.25 if i == j else polarized[i] * polarized[j] for j in range(6)] for i in range(6)]
    for i, j in ((2, 4), (3, 5)):
        zz[i][j] = zz[j][i] = -0.25
    connected = _connected_zz_windows(geometry, polarized, {"correlations": {"zz_real": zz}})
    central_right = next(window["statistics"] for window in connected["windows"]
                         if window["left_window"] == "central" and window["right_window"] == "right_edge")
    assert central_right["count"] == 4 and central_right["mean"] == -0.125
    assert math.isclose(central_right["rms"], 1 / (4 * math.sqrt(2)), abs_tol=2e-15)
    assert connected["fixed_Q_connected_row_sum_max_abs"] == 0
    assert _connected_zz_windows(geometry, polarized, {})["status"] == "not_recorded"
    assert _connected_zz_windows(geometry, polarized, {"correlations": {"zz_real": [[0.0]]}})["status"] == "invalid_recorded_matrix"

    lower = [1 / 54] * 27
    column_delta = [-0.2, 0.5, 0.7]
    higher = [value + column_delta[site["x"]] / 9 for value, site in zip(lower, fourier_geometry["sites"])]
    delta = _signed_column_delta(fourier_geometry, lower, higher, 2)
    assert all(math.isclose(row["signed_delta_sz"], expected, abs_tol=2e-15) for row, expected in zip(delta["columns"], column_delta))
    assert delta["columns"][0]["signed_delta_sz"] < 0 and delta["sum_rule_passed"]
    highest = [value + change / 9 for value, change in zip(higher, [-0.1] * 9 + [0.2] * 9 + [0.9] * 9)]
    sector_cases = []
    for q, values in ((1, lower), (3, higher), (5, highest)):
        case = artificial_case(q, -10 + q / 10)
        case["configuration"] = {**fourier_geometry, "Q": q}
        case["solver"], case["backend_source_sha256"] = {"maxdim": 256}, {}
        case["batches"][0]["sz_profile"] = values
        sector_cases.append(case)
    sector_changes = _sector_column_changes(sector_cases)
    assert len(sector_changes["comparisons"]) == 2
    assert all(comparison["delta"]["sum_rule_passed"] for comparison in sector_changes["comparisons"])
    # Cover nested array-of-table serialization used by cases/batches/columns.
    payload = {"passed": True, "value": 0.25, "rows": [{"a": [1, 2], "inner": {"x": "λ"}},
                                                       {"a": [3], "children": [{"b": False}]}]}
    assert tomllib.loads(dumps_toml(payload)) == payload
    return {"passed": True, "elapsed_seconds": time.monotonic() - start,
            "checks": ["RMS formulas", "incomplete/unverified batch exclusion", "Q=2M field differences", "candidate spread range", "missing-sector handling",
                       "execution allowlist and companion filenames", "audit digest linkage and mismatch", "strict nonempty recorded audit checks",
                       "known circumference site/bond shift and one-to-one matching",
                       "stationarity remains separate", "coordinate and endpoint weighting", "one-sweep missing measurement preservation",
                       "chi-change stationarity exclusion", "legacy comparison scope", "HdaggerH imaginary coverage",
                       "hand Fourier coefficients and phase sign", "connected Czz window arithmetic and missing data",
                       "signed Q1->3 and Q3->5 column sums", "nested TOML round-trip"],
            "numerical_solver_run": False}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", nargs="?", type=Path)
    parser.add_argument("inputs", nargs="*", type=Path)
    parser.add_argument("--selfcheck", action="store_true")
    args = parser.parse_args()
    if args.selfcheck:
        print(json.dumps(selfcheck(), ensure_ascii=False))
        return 0
    if args.output is None or not args.inputs:
        parser.error("provide OUTPUT.toml and one or more validation inputs")
    if args.output.suffix.lower() != ".toml":
        parser.error("output must be a .toml analysis artifact")
    source_paths = {source.resolve() for path in args.inputs for source in (path, *_companion_paths(path))}
    if args.output.resolve() in source_paths:
        parser.error("output must not overwrite an input validation or companion provenance file")
    report = summarize(args.inputs)
    content = dumps_toml(report)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    temporary = None
    try:
        with tempfile.NamedTemporaryFile(mode="w", encoding="utf-8", dir=args.output.parent,
                                         prefix=f".{args.output.name}.", suffix=".tmp", delete=False) as stream:
            temporary = stream.name
            stream.write(content)
        os.replace(temporary, args.output)
    finally:
        if temporary and os.path.exists(temporary):
            os.unlink(temporary)
    print(json.dumps({"output": str(args.output.resolve()), "status": report["status"],
                      "eligible_cases": report["eligible_case_count"], "eligible_batches": report["eligible_batch_count"]}))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
