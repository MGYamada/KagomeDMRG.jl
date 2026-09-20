#!/usr/bin/env python3
"""Analyze saved NN N18 Q2 ED vectors within an inclusive 120 second limit."""

import argparse
import math
import os
from pathlib import Path
import sys

sys.dont_write_bytecode = True
from run_static27 import supervise_with_resources


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--wall-seconds", type=float, default=120.0)
    args = parser.parse_args()
    if not math.isfinite(args.wall_seconds) or not 0 < args.wall_seconds <= 120:
        parser.error("--wall-seconds must be finite, positive, and at most 120")
    if os.name != "posix" or not (sys.platform == "darwin" or sys.platform.startswith("linux")):
        parser.error("this launcher requires Darwin or Linux resource accounting")
    command = ["julia", "--project=research", "--startup-file=no", "--threads=1",
               "examples/analyze_static18_neutral.jl", str(args.output.resolve())]
    environment = os.environ.copy()
    for key in ("JULIA_NUM_THREADS", "JULIA_NUM_GC_THREADS", "JULIA_NUM_PRECOMPILE_TASKS",
                "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS"):
        environment[key] = "1"
    try:
        return supervise_with_resources(command, args.output, args.wall_seconds, env=environment)
    except (OSError, ValueError) as error:
        print(f"static18 neutral launcher: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
