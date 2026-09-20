#!/usr/bin/env python3
"""Measure the saved CSL N36 trial's missing moment within 300 seconds."""

import argparse
import math
import os
from pathlib import Path
import sys

sys.dont_write_bytecode = True
from run_static27 import supervise_with_resources

ROOT = Path(__file__).resolve().parent.parent
BACKEND = ROOT / "outputs/source-snapshots/csl36-pre-env-split-20260920"


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--wall-seconds", type=float, default=300.0)
    args = parser.parse_args()
    if not math.isfinite(args.wall_seconds) or not 0 < args.wall_seconds <= 300:
        parser.error("--wall-seconds must be finite, positive, and at most 300")
    if os.name != "posix" or not (sys.platform == "darwin" or sys.platform.startswith("linux")):
        parser.error("requires Darwin or Linux process/resource accounting")
    if not (BACKEND / "Project.toml").is_file():
        parser.error("the archived backend project is required")
    command = ["julia", "--project=.", "--startup-file=no", "--threads=1",
               str(ROOT / "examples/audit_csl_saved_moment.jl"), str(args.output.resolve())]
    environment = os.environ.copy()
    for key in ("JULIA_NUM_THREADS", "JULIA_NUM_GC_THREADS", "JULIA_NUM_PRECOMPILE_TASKS",
                "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS"):
        environment[key] = "1"
    try:
        return supervise_with_resources(command, args.output, args.wall_seconds,
                                        cwd=BACKEND, env=environment)
    except (OSError, ValueError) as error:
        print(f"CSL saved moment launcher: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
