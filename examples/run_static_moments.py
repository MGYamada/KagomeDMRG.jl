#!/usr/bin/env python3
"""Audit saved N27 moments using the original backend within 120 seconds."""

import argparse
import math
import os
from pathlib import Path
import sys

sys.dont_write_bytecode = True
from run_static27 import supervise_with_resources

ROOT = Path(__file__).resolve().parent.parent
BACKEND = ROOT / "outputs/source-snapshots/n27-refine-before-chirality-20260919"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path, nargs="?",
                        default=ROOT / "outputs/p4-campaign-moment-audit")
    parser.add_argument("--wall-seconds", type=float, default=120.0)
    args = parser.parse_args()
    if not math.isfinite(args.wall_seconds) or not 0 < args.wall_seconds <= 120:
        parser.error("--wall-seconds must be finite, positive, and at most 120")
    if os.name != "posix" or not (sys.platform == "darwin" or sys.platform.startswith("linux")):
        parser.error("requires Darwin or Linux process/resource accounting")
    if not (BACKEND / "Project.toml").is_file():
        parser.error("the archived backend project is required")
    command = ["julia", "--project=.", "--startup-file=no", "--threads=1",
               str(ROOT / "examples/audit_static_moments.jl"), str(args.output.resolve())]
    environment = os.environ.copy()
    for key in ("JULIA_NUM_THREADS", "JULIA_NUM_GC_THREADS", "JULIA_NUM_PRECOMPILE_TASKS",
                "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS"):
        environment[key] = "1"
    try:
        return supervise_with_resources(command, args.output, args.wall_seconds,
                                        cwd=BACKEND, env=environment)
    except (OSError, ValueError) as error:
        print(f"static moments launcher: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
