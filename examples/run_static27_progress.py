#!/usr/bin/env python3
"""Run the N27 two-sweep progress study within one inclusive 600 s budget."""

import argparse
import math
import os
from pathlib import Path
import sys

sys.dont_write_bytecode = True
from run_static27 import MAX_WALL_SECONDS, supervise_with_resources


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--wall-seconds", type=float, default=MAX_WALL_SECONDS)
    args = parser.parse_args()
    if not math.isfinite(args.wall_seconds) or not 0 < args.wall_seconds <= MAX_WALL_SECONDS:
        parser.error("--wall-seconds must be finite, positive, and at most 600")
    if os.name != "posix" or not (sys.platform == "darwin" or sys.platform.startswith("linux")):
        parser.error("this launcher requires Darwin or Linux process/resource accounting")
    output = args.output.resolve()
    command = ["julia", "--project=research", "--startup-file=no", "--threads=1",
               "examples/validate_static27_progress.jl", str(output)]
    environment = os.environ.copy()
    for key in ("JULIA_NUM_THREADS", "JULIA_NUM_GC_THREADS", "JULIA_NUM_PRECOMPILE_TASKS",
                "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS"):
        environment[key] = "1"
    try:
        return supervise_with_resources(command, output, args.wall_seconds, env=environment)
    except (OSError, ValueError) as error:
        print(f"static27 progress launcher: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
