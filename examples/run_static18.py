#!/usr/bin/env python3
"""Run the static N=18 worker with an external, inclusive wall-time limit.

Usage: python3 examples/run_static18.py OUTPUT [--wall-seconds 600]
The final <=2 seconds of the limit are reserved for process-group termination.
Only execution.toml is owned here; a worker's validation.toml is never promoted
to success merely because its process exited successfully.
"""

import argparse
import json
import math
import os
from pathlib import Path
import signal
import subprocess
import sys
import tempfile
import time


ROOT = Path(__file__).resolve().parent.parent
MAX_WALL_SECONDS = 600.0


def write_record(path, record):
    """Atomically replace a small flat TOML allowlist, using only stdlib."""
    def encode(value):
        if isinstance(value, bool):
            return "true" if value else "false"
        if isinstance(value, str):
            # Keep Unicode scalars intact; TOML also requires DEL to be escaped.
            return json.dumps(value, ensure_ascii=False).replace("\x7f", "\\u007f")
        if isinstance(value, list):
            return "[" + ", ".join(encode(item) for item in value) + "]"
        return str(value)

    temporary = None
    try:
        with tempfile.NamedTemporaryFile(
            mode="w", encoding="utf-8", dir=path.parent,
            prefix=".execution-", suffix=".toml", delete=False,
        ) as stream:
            temporary = Path(stream.name)
            for key, value in record.items():
                stream.write(f"{key} = {encode(value)}\n")
            stream.flush()
            os.fsync(stream.fileno())
        os.replace(temporary, path)
    finally:
        if temporary is not None and temporary.exists():
            temporary.unlink()


def signal_group(process, signum):
    # start_new_session=True below makes this process the group leader. Never
    # use name-based process lookup or signal any pre-existing Julia process.
    try:
        os.killpg(process.pid, signum)
    except ProcessLookupError:
        pass


def supervise(command, output, wall_seconds, *, cwd=ROOT, env=None):
    """Start exactly one worker; private function also supports tiny mock runs."""
    output = Path(output).resolve()
    output.mkdir(parents=True, exist_ok=True)
    if any(output.iterdir()):
        raise ValueError("output directory must be new or empty")
    execution = output / "execution.toml"
    # Exclusive creation also rejects another launcher claiming the same dir.
    with execution.open("x", encoding="utf-8"):
        pass
    grace = min(2.0, wall_seconds / 5.0)
    # Reserve a small part of the grace for reaping after SIGKILL, so even a
    # worker ignoring SIGTERM is killed before the inclusive deadline.
    reap_reserve = min(0.1, grace / 4.0)
    record = {
        "schema_version": 1,
        "status": "starting",
        "command": command,
        "wall_limit_seconds": wall_seconds,
        "termination_grace_seconds": grace,
        "termination_grace_included_in_wall_limit": True,
        "julia_threads": 1,
        "blas_threads": 1,
        "worker_exit_confirmed": False,
        "worker_validation_status": "not_interpreted",
    }
    write_record(execution, record)
    process = None
    # Measured immediately before Popen: startup and all Julia compilation are
    # charged to the same budget as the solves.
    started = time.monotonic()
    deadline = started + wall_seconds
    try:
        process = subprocess.Popen(
            command, cwd=cwd, env=env, start_new_session=True,
        )
        record["status"] = "running"
        write_record(execution, record)
        try:
            process.wait(timeout=max(0.0, deadline - grace - time.monotonic()))
            record["status"] = "exited" if process.returncode == 0 else "failed"
        except subprocess.TimeoutExpired:
            record["status"] = "timed_out"
            record["sigterm_elapsed_seconds"] = time.monotonic() - started
            signal_group(process, signal.SIGTERM)
            try:
                process.wait(timeout=max(0.0, deadline - reap_reserve - time.monotonic()))
            except subprocess.TimeoutExpired:
                pass
            # Always clean up our group after a timeout: the direct child may
            # have exited on SIGTERM while one of its descendants did not.
            signal_group(process, signal.SIGKILL)
            record["sigkill_elapsed_seconds"] = time.monotonic() - started
            try:
                process.wait(timeout=max(0.0, deadline - time.monotonic()))
            except subprocess.TimeoutExpired:
                pass
    except KeyboardInterrupt:
        record["status"] = "interrupted"
    except OSError as error:
        record["status"] = "launch_or_io_failed"
        record["error_errno"] = error.errno or 0
    finally:
        if process is not None:
            if process.poll() is None:
                signal_group(process, signal.SIGKILL)
                try:
                    process.wait(timeout=max(0.0, deadline - time.monotonic()))
                except subprocess.TimeoutExpired:
                    pass
            record["worker_exit_confirmed"] = process.poll() is not None
            if process.returncode is not None:
                record["worker_exit_code"] = process.returncode
        record["wall_elapsed_seconds"] = time.monotonic() - started
        write_record(execution, record)
    print(f"static18 launcher: {record['status']} after "
          f"{record['wall_elapsed_seconds']:.3f} s", flush=True)
    if record["status"] == "timed_out":
        return 124
    if record["status"] == "interrupted":
        return 130
    return 0 if record["status"] == "exited" else 1


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    parser.add_argument("--wall-seconds", type=float, default=MAX_WALL_SECONDS,
                        help="inclusive runtime limit, >0 and <=600 (default 600)")
    args = parser.parse_args()
    if not math.isfinite(args.wall_seconds) or not 0 < args.wall_seconds <= MAX_WALL_SECONDS:
        parser.error("--wall-seconds must be finite, positive, and at most 600")
    if os.name != "posix":
        parser.error("this launcher requires POSIX process groups")
    output = args.output.resolve()
    command = ["julia", "--project=research", "--startup-file=no", "--threads=1",
               "examples/validate_static18.jl", str(output)]
    environment = os.environ.copy()
    for key in ("JULIA_NUM_THREADS", "JULIA_NUM_GC_THREADS", "JULIA_NUM_PRECOMPILE_TASKS",
                "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS", "OMP_NUM_THREADS",
                "VECLIB_MAXIMUM_THREADS", "BLIS_NUM_THREADS"):
        environment[key] = "1"
    try:
        return supervise(command, output, args.wall_seconds, env=environment)
    except (OSError, ValueError) as error:
        # No environment or connection diagnostics are persisted.
        print(f"static18 launcher: {error}", file=sys.stderr)
        return 2


if __name__ == "__main__":
    sys.exit(main())
