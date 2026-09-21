#!/usr/bin/env python3
"""Run one complete, explicitly budgeted preparation comparison, one worker/thread."""
import argparse
import os
from pathlib import Path
import sys
import tomllib

sys.dont_write_bytecode = True
from run_static27 import supervise_with_resources

ROOT = Path(__file__).resolve().parent.parent


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--config", type=Path,
                        default=ROOT / "examples/configs/vbc54_comparison.toml")
    args = parser.parse_args()
    cfg = tomllib.loads(args.config.read_text())
    # This is a budget for all three branches, not a per-sweep timeout.
    if cfg["wall_seconds"] != 2700 or cfg["branches"] != ["random", "hourglass", "windmill"]:
        parser.error("expected the declared three-branch, 2700 second comparison")
    env = os.environ.copy()
    for key in ("JULIA_NUM_THREADS", "JULIA_NUM_GC_THREADS", "JULIA_NUM_PRECOMPILE_TASKS",
                "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS"):
        env[key] = "1"
    command = ["julia", "--project=research", "--startup-file=no", "--threads=1",
               str(ROOT / "examples/research_vbc54.jl"), str(args.output.resolve()),
               str(args.config.resolve())]
    return supervise_with_resources(command, args.output, cfg["wall_seconds"], cwd=ROOT, env=env)


if __name__ == "__main__":
    sys.exit(main())
