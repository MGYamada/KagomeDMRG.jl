#!/usr/bin/env python3
"""Run one NN sweep (<=900 s) or its saved-state diagnostics (<=180 s).

Each stage requires a new output directory. The original solve record is never
rewritten by diagnostics. Run stages sequentially, with one numerical worker.
"""

import argparse
import hashlib
import math
import os
from pathlib import Path
import sys
import tomllib

sys.dont_write_bytecode = True
from run_static27 import supervise_with_resources

ROOT = Path(__file__).resolve().parent.parent
BACKEND = ROOT / "outputs/source-snapshots/n27-refine-before-chirality-20260919"
LIMITS = {"solve": 900.0, "diagnose": 180.0}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("stage", choices=tuple(LIMITS))
    parser.add_argument("config", type=Path)
    parser.add_argument("output", type=Path)
    parser.add_argument("--backend", type=Path, default=BACKEND)
    parser.add_argument("--solve-output", type=Path)
    parser.add_argument("--wall-seconds", type=float)
    args = parser.parse_args()
    limit = LIMITS[args.stage]
    wall = limit if args.wall_seconds is None else args.wall_seconds
    if not math.isfinite(wall) or not 0 < wall <= limit:
        parser.error(f"--wall-seconds must be finite, positive, and at most {limit:g}")
    if os.name != "posix" or not (sys.platform == "darwin" or sys.platform.startswith("linux")):
        parser.error("requires Darwin or Linux process/resource accounting")
    backend, config, output = args.backend.resolve(), args.config.resolve(), args.output.resolve()
    if backend != BACKEND.resolve() or not (backend / "Project.toml").is_file():
        parser.error("the original archived N27 backend is required")
    if not config.is_file():
        parser.error("config must exist")
    if (args.stage == "diagnose") != (args.solve_output is not None):
        parser.error("--solve-output is required only for diagnose")
    command = ["julia", "--project=.", "--startup-file=no", "--threads=1",
               str(ROOT / "examples/research_static_split.jl"), args.stage,
               str(output), str(config)]
    if args.stage == "diagnose":
        solve = args.solve_output.resolve()
        try:
            validation_bytes = (solve / "validation.toml").read_bytes()
            execution_bytes = (solve / "execution.toml").read_bytes()
            validation = tomllib.loads(validation_bytes.decode("utf-8"))
            execution = tomllib.loads(execution_bytes.decode("utf-8"))
            if validation.get("status") != "completed_solve_diagnostics_deferred":
                raise ValueError("a completed split solve is required")
            if (execution.get("status") != "exited" or
                    execution.get("worker_exit_code") != 0 or
                    execution.get("worker_exit_confirmed") is not True):
                raise ValueError("successful solve worker exit is required")
        except (OSError, ValueError) as error:
            parser.error(str(error))
        command.extend([str(solve), hashlib.sha256(validation_bytes).hexdigest(),
                        hashlib.sha256(execution_bytes).hexdigest()])
    environment = os.environ.copy()
    for key in ("JULIA_NUM_THREADS", "JULIA_NUM_GC_THREADS", "JULIA_NUM_PRECOMPILE_TASKS",
                "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS"):
        environment[key] = "1"
    try:
        # The shared supervisor accepts an explicit inclusive bound. Its older
        # <=600 restriction belongs to the old CLI, which is not invoked here.
        return supervise_with_resources(command, output, wall, cwd=backend, env=environment)
    except (OSError, ValueError) as error:
        print(f"static split launcher: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
