#!/usr/bin/env python3
"""Run one explicit static/CSL case, one thread, with a <=600 s inclusive limit."""
import argparse
import math
import os
from pathlib import Path
import sys

sys.dont_write_bytecode = True
from run_static27 import supervise_with_resources
ROOT = Path(__file__).resolve().parent.parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("worker", choices=("static", "csl"))
    parser.add_argument("config", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--backend", type=Path, default=ROOT)
    parser.add_argument("--wall-seconds", type=float, default=600.0)
    args = parser.parse_args()
    if not math.isfinite(args.wall_seconds) or not 0 < args.wall_seconds <= 600:
        parser.error("--wall-seconds must be finite, positive, and at most 600")
    if os.name != "posix" or not (sys.platform == "darwin" or sys.platform.startswith("linux")):
        parser.error("requires Darwin or Linux process/resource accounting")
    backend, config, output = args.backend.resolve(), args.config.resolve(), args.output.resolve()
    if not (backend / "Project.toml").is_file() or not config.is_file():
        parser.error("backend project and config file must exist")
    # Current checkouts keep their pinned environment in research; archived
    # backends without that layout retain their original root environment.
    project = "research" if (backend / "research" / "Project.toml").is_file() else "."
    command = ["julia", f"--project={project}", "--startup-file=no", "--threads=1",
               str(ROOT / "examples" / f"research_{args.worker}.jl"), str(output), str(config)]
    environment = os.environ.copy()
    for key in ("JULIA_NUM_THREADS", "JULIA_NUM_GC_THREADS", "JULIA_NUM_PRECOMPILE_TASKS",
                "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS"):
        environment[key] = "1"
    try:
        return supervise_with_resources(command, output, args.wall_seconds, cwd=backend, env=environment)
    except (OSError, ValueError) as error:
        print(f"research launcher: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
