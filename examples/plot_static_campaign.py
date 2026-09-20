#!/usr/bin/env python3
"""Plot the latest completed N27/N54 density trials without translating profiles.

Uses three immutable validation byte snapshots. Requires Matplotlib; no Julia,
DMRG, or checkpoint deserialization. The whole script has a 120 second limit.
"""

from __future__ import annotations

import argparse
import datetime as dt
import hashlib
import json
import math
import os
from pathlib import Path
import signal
import sys
import tempfile
import time
import tomllib

ROOT = Path(__file__).resolve().parent.parent
DEFAULT_INPUTS = [ROOT / "outputs" / case / "validation.toml" for case in (
    "p4-campaign-nn27-period9-chi128", "p4-campaign-nn27-period27-chi128",
    "p4-campaign-nn54-period27-chi128")]
DEFAULT_STEM = ROOT / "docs/research/figures/p4_campaign_density"


def sha256(payload: bytes) -> str:
    return hashlib.sha256(payload).hexdigest()


def select_trial(path: Path, expected_n: int, expected_kind: str) -> dict:
    payload = path.read_bytes()
    record = tomllib.loads(payload.decode("utf-8"))
    config, geometry = record["config"], record["configuration"]
    if (record["N"], config["maxdim"], config["initialization"]) != (expected_n, 128, expected_kind):
        raise ValueError("expected the declared N27/N54 chi128 ordered-seed cases")
    if (geometry["Lx"], geometry["Ly"], record["Q"]) != (expected_n // 9, 3, expected_n // 9):
        raise ValueError("expected Ly=3 and Q=N/9")
    if record["theta"] != 0 or any(geometry["hz"]):
        raise ValueError("expected zero-flux, zero-field profiles")
    if any(record.get(key) is False for key in ("sources_unchanged", "backend_unchanged", "parent_unchanged")):
        raise ValueError("input records an integrity failure")
    if any(b["family"] != "J1" or b["Jxy"] != 1 or b["Jz"] != 1 for b in geometry["bonds"]):
        raise ValueError("expected isotropic nearest-neighbor J=1")
    rows = [row for row in record["batches"]
            if row.get("status") == "diagnostics_passed_accuracy_separate"
            and row.get("active_phase") == "completed"
            and row.get("integrity_passed") is True
            and row.get("checkpoint_verified") is True]
    if not rows:
        raise ValueError(f"no completed integrity-verified trial in {path.name}")
    selected = max(rows, key=lambda row: row["cumulative_sweeps"])
    profile = selected["sz_profile"]
    lx, ly = geometry["Lx"], geometry["Ly"]
    if len(profile) != expected_n or not all(math.isfinite(v) and abs(v) <= 0.5 + 1e-10 for v in profile):
        raise ValueError("invalid physical Sz profile")
    if len(geometry["sites"]) != expected_n:
        raise ValueError("geometry site count disagrees with N")
    for index, site in enumerate(geometry["sites"]):
        x, remainder = divmod(index, 3 * ly)
        y, sub = divmod(remainder, 3)
        if (site["index"], site["x"], site["y"], site["sublattice"]) != (index + 1, x, y, "ABC"[sub]):
            raise ValueError("unsupported site ordering")
    columns = [math.fsum(profile[9 * x:9 * (x + 1)]) for x in range(lx)]
    if abs(math.fsum(columns) - record["Q"] / 2) > 1e-10:
        raise ValueError("profile fails total-charge sum rule")
    if len(selected["column_sz"]) != lx or max(abs(a - b) for a, b in zip(columns, selected["column_sz"])) > 1e-10:
        raise ValueError("stored column sums disagree with the site profile")
    variance = selected["variance"]
    if not math.isfinite(variance) or variance < -selected["variance_roundoff_scale"]:
        raise ValueError("invalid recorded variance")
    return {
        "input_path": os.path.relpath(path.resolve(), ROOT), "input_sha256": sha256(payload),
        "input_bytes": len(payload), "case_id": record["case_id"],
        "worker_status_recorded": record["status"],
        "selected_batch": selected["batch"], "selected_cumulative_sweeps": selected["cumulative_sweeps"],
        "selection": "latest completed integrity-verified batch in this byte snapshot",
        "N": expected_n, "Q": record["Q"], "Lx": lx, "Ly": ly,
        "initialization": expected_kind, "configured_maxdim": config["maxdim"],
        "actual_maxlinkdim": selected["maxlinkdim"], "variance": variance,
        "variance_per_site": variance / expected_n, "energy": selected["energy"],
        "checkpoint": selected["checkpoint"],
        "checkpoint_metadata_sha256_recorded": selected["checkpoint_metadata_sha256"],
        "checkpoint_payload_sha256_recorded": selected["checkpoint_payload_sha256"],
        "sz_profile": profile, "column_sz": columns,
        "uniform_column_reference": record["Q"] / (2 * lx),
        "row_order": [f"{y}/{sub}" for y in range(ly) for sub in "ABC"],
        "profile_translation_applied": False,
    }


def render(trials: list[dict], stem: Path) -> dict:
    # Set these before importing Matplotlib/NumPy. No persistent user cache.
    for key in ("OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS", "NUMEXPR_NUM_THREADS"):
        os.environ[key] = "1"
    with tempfile.TemporaryDirectory(prefix="kagome-campaign-mpl-", dir="/tmp") as cache:
        os.environ["MPLCONFIGDIR"] = cache
        os.environ["XDG_CACHE_HOME"] = cache
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.colors import Normalize

        plt.rcParams.update({"font.family": "DejaVu Sans", "font.size": 10,
                             "axes.titlesize": 11, "axes.labelsize": 10,
                             "svg.fonttype": "none", "savefig.facecolor": "white"})
        fig = plt.figure(figsize=(12.5, 7.1), facecolor="white")
        grid = fig.add_gridspec(2, 4, width_ratios=[3, 3, 6, 0.28], height_ratios=[9, 3.3],
                               left=0.085, right=0.92, bottom=0.19, top=0.73,
                               wspace=0.24, hspace=0.30)
        norm = Normalize(vmin=-0.5, vmax=0.5)
        all_columns = [value for trial in trials for value in trial["column_sz"]] + [0.5]
        low, high = min(all_columns), max(all_columns)
        padding = max(0.025, 0.15 * (high - low))
        heatmaps = []
        for index, trial in enumerate(trials):
            lx = trial["Lx"]
            ax = fig.add_subplot(grid[0, index])
            matrix = [[trial["sz_profile"][9 * x + r] for x in range(lx)] for r in range(9)]
            heatmap = ax.imshow(matrix, cmap="RdBu_r", norm=norm, interpolation="nearest",
                                origin="upper", aspect="auto", extent=(-0.5, lx - 0.5, 8.5, -0.5))
            heatmaps.append(heatmap)
            ax.set_xticks(range(lx))
            ax.tick_params(axis="x", labelbottom=False, bottom=False)
            ax.set_yticks(range(9), trial["row_order"] if index == 0 else [""] * 9)
            ax.tick_params(axis="y", length=0)
            if index == 0:
                ax.set_ylabel("(y, sublattice)", labelpad=10)
            ax.set_xticks([x - 0.5 for x in range(lx + 1)], minor=True)
            ax.set_yticks([r - 0.5 for r in range(10)], minor=True)
            ax.grid(which="minor", color="white", linewidth=0.55, alpha=0.75)
            ax.tick_params(which="minor", bottom=False, left=False)
            ax.set_title(f"N = {trial['N']} | {trial['initialization']} seed\n"
                         f"χ = {trial['configured_maxdim']} | {trial['selected_cumulative_sweeps']} sweeps\n"
                         f"Var(H) = {trial['variance']:.4g} J²", pad=11, linespacing=1.4)

            col = fig.add_subplot(grid[1, index])
            col.axhline(trial["uniform_column_reference"], color="0.5", linestyle="--", linewidth=1.0)
            col.plot(range(lx), trial["column_sz"], color="0.12", marker="o", markersize=4,
                     linewidth=1.25)
            col.set_xlim(-0.5, lx - 0.5)
            col.set_ylim(low - padding, high + padding)
            col.set_xticks(range(lx))
            col.set_xlabel("open-axis column x")
            col.spines[["top", "right"]].set_visible(False)
            col.grid(axis="y", color="0.9", linewidth=0.6)
            if index == 0:
                col.set_ylabel("Column Σ Sz")
            else:
                col.tick_params(axis="y", labelleft=False)

        colorbar = fig.colorbar(heatmaps[0], cax=fig.add_subplot(grid[0, 3]), ticks=[-0.5, 0, 0.5])
        colorbar.set_label("Physical ⟨Sz⟩", labelpad=10)
        fig.suptitle("Nearest-neighbor 1/9 density profiles", x=0.085, y=0.98,
                     ha="left", fontsize=15, fontweight="normal")
        fig.text(0.085, 0.935, "Finite trials; not converged. Latest completed batch from each saved input.",
                 ha="left", fontsize=10.5, color="0.25")
        fig.text(0.085, 0.090, "Common color scale: −0.5 to +0.5 physical Sz. Raw coordinate profiles; no translations.",
                 ha="left", fontsize=10)
        fig.text(0.085, 0.057, "Dashed line: uniform 1/9 reference = 0.5 per column. N54 adjacent sectors are not measured.",
                 ha="left", fontsize=10)
        stem.parent.mkdir(parents=True, exist_ok=True)
        png, svg = stem.with_suffix(".png"), stem.with_suffix(".svg")
        fig.savefig(png, dpi=180, metadata={"Software": "KagomeDMRG.jl plot_static_campaign.py"})
        fig.savefig(svg, metadata={"Title": "Nearest-neighbor 1/9 finite-trial density profiles",
                                  "Creator": "KagomeDMRG.jl plot_static_campaign.py", "Date": None})
        plt.close(fig)
        return {"matplotlib_version": matplotlib.__version__,
                "png": {"path": os.path.relpath(png, ROOT), "sha256": sha256(png.read_bytes())},
                "svg": {"path": os.path.relpath(svg, ROOT), "sha256": sha256(svg.read_bytes())}}


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--inputs", type=Path, nargs=3, default=DEFAULT_INPUTS,
                        metavar=("N27_PERIOD9", "N27_PERIOD27", "N54_PERIOD27"))
    parser.add_argument("--output-stem", type=Path, default=DEFAULT_STEM)
    args = parser.parse_args()
    start = time.monotonic()
    def timeout(_signum, _frame):
        raise TimeoutError("plotting exceeded the 120 second wall limit")
    signal.signal(signal.SIGALRM, timeout)
    signal.setitimer(signal.ITIMER_REAL, 120.0)
    try:
        trials = [select_trial(path, n, kind) for path, (n, kind) in
                  zip(args.inputs, ((27, "period9"), (27, "period27"), (54, "period27")))]
        stem = args.output_stem.resolve()
        artifacts = render(trials, stem)
        provenance = {"schema_version": 1, "recorded_at_utc": dt.datetime.now(dt.timezone.utc).isoformat(),
                      "script": "examples/plot_static_campaign.py", "script_sha256": sha256(Path(__file__).read_bytes()),
                      "python_version": sys.version.split()[0], "wall_limit_seconds": 120,
                      "wall_elapsed_seconds": time.monotonic() - start, "new_numerical_solver_run": False,
                      "color_scale": {"minimum": -0.5, "maximum": 0.5, "units": "physical Sz", "shared": True},
                      "no_coordinate_shift": True, "trial_status": "finite_trials_not_converged",
                      "N54_adjacent_sectors": "not_measured_in_this_comparison", "inputs": trials,
                      "artifacts": artifacts}
        stem.with_suffix(".provenance.json").write_text(json.dumps(provenance, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")
        print(json.dumps({"png": artifacts["png"]["path"], "svg": artifacts["svg"]["path"],
                          "selected_sweeps": [t["selected_cumulative_sweeps"] for t in trials],
                          "wall_seconds": provenance["wall_elapsed_seconds"]}))
    finally:
        signal.setitimer(signal.ITIMER_REAL, 0.0)
    return 0


if __name__ == "__main__":
    sys.exit(main())
