#!/usr/bin/env python3
"""Run the single static N=27 worker with a <=600 s inclusive wall limit.

Usage: python3 examples/run_static27.py OUTPUT [--wall-seconds 600]
Requires Python >=3.11 on Darwin or Linux. The N=18 supervisor supplies the
same process-group timeout and empty-output protection without modification.
Only execution.toml is updated here; worker validation claims are untouched.
"""

import argparse
import math
import os
from pathlib import Path
import resource
import sys
import tomllib

# Import the existing supervisor without adding cache artifacts to the examples.
sys.dont_write_bytecode = True
from run_static18 import MAX_WALL_SECONDS, ROOT, supervise, write_record


def resource_record(before, after, worker_exit_confirmed, platform=sys.platform):
    """Return explicit child accounting; never subtract two RSS high waters."""
    if platform == "darwin":
        unit, multiplier = "bytes", 1
    elif platform.startswith("linux"):
        unit, multiplier = "KiB", 1024
    else:
        raise ValueError("ru_maxrss conversion is supported only on Darwin and Linux")
    prior_child_usage = any(
        getattr(before, key) != 0 for key in ("ru_maxrss", "ru_utime", "ru_stime")
    )
    record = {
        "resource_measurement_method": "resource.getrusage(resource.RUSAGE_CHILDREN)",
        "resource_measurement_scope": "single worker reaped by this launcher; includes descendant usage only as accounted by the OS",
        "resource_excludes_launcher": True,
        "resource_is_simultaneous_process_tree_peak": False,
        "resource_prior_child_usage_detected": prior_child_usage,
        "resource_cpu_method": "difference of cumulative child ru_utime and ru_stime before launch and after exit",
        "resource_peak_rss_method": "post-exit ru_maxrss high water; available only with no prior child usage",
        "resource_ru_maxrss_native_unit": unit,
        "resource_ru_maxrss_bytes_per_native_unit": multiplier,
    }
    if not worker_exit_confirmed:
        record["resource_measurement_status"] = "unavailable_worker_exit_not_confirmed"
        return record
    user_seconds = after.ru_utime - before.ru_utime
    system_seconds = after.ru_stime - before.ru_stime
    record.update({
        "worker_cpu_user_seconds": user_seconds,
        "worker_cpu_system_seconds": system_seconds,
        "worker_cpu_total_seconds": user_seconds + system_seconds,
    })
    if prior_child_usage:
        # RUSAGE_CHILDREN is process-lifetime accounting. A previous child's
        # peak cannot be removed from a maximum to recover this worker's peak.
        record["resource_measurement_status"] = "cpu_captured_peak_rss_unavailable_prior_child_usage"
    else:
        record["resource_measurement_status"] = "captured"
        record["worker_peak_rss_bytes"] = int(after.ru_maxrss * multiplier)
    return record


def supervise_with_resources(command, output, wall_seconds, *, cwd=ROOT, env=None):
    before = resource.getrusage(resource.RUSAGE_CHILDREN)
    result = supervise(command, output, wall_seconds, cwd=cwd, env=env)
    after = resource.getrusage(resource.RUSAGE_CHILDREN)
    execution = Path(output).resolve() / "execution.toml"
    with execution.open("rb") as stream:
        record = tomllib.load(stream)
    record.update(resource_record(before, after, record["worker_exit_confirmed"]))
    write_record(execution, record)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--wall-seconds", type=float, default=MAX_WALL_SECONDS,
                        help="inclusive runtime limit, >0 and <=600 (default 600)")
    args = parser.parse_args()
    if not math.isfinite(args.wall_seconds) or not 0 < args.wall_seconds <= MAX_WALL_SECONDS:
        parser.error("--wall-seconds must be finite, positive, and at most 600")
    if os.name != "posix" or not (sys.platform == "darwin" or sys.platform.startswith("linux")):
        parser.error("this launcher requires Darwin or Linux process/resource accounting")
    output = args.output.resolve()
    command = ["julia", "--project=research", "--startup-file=no", "--threads=1",
               "examples/validate_static27.jl", str(output)]
    environment = os.environ.copy()
    for key in ("JULIA_NUM_THREADS", "JULIA_NUM_GC_THREADS", "JULIA_NUM_PRECOMPILE_TASKS",
                "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS"):
        environment[key] = "1"
    try:
        return supervise_with_resources(command, output, args.wall_seconds, env=environment)
    except (OSError, ValueError) as error:
        print(f"static27 launcher: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
