# SPDX-License-Identifier: MIT
"""
Connect the container's cluster-generic executor to host Slurm commands.

The Sophia image need not contain sbatch/sacct/scancel. Requests and replies use
the shared project filesystem; the host build allocation runs the small bridge.
"""

import os
import sys
import time
import uuid
from pathlib import Path

from snakemake.utils import read_job_properties


def request(operation, arguments):
    spool = Path(os.environ["BORNHOLM_WORKQ_BRIDGE"])
    token = spool / uuid.uuid4().hex
    payload = [operation, *map(str, arguments)]
    token.with_suffix(".tmp").write_bytes(("\0".join(payload) + "\0").encode())
    token.with_suffix(".tmp").replace(token.with_suffix(".request"))
    deadline = time.monotonic() + 120
    while not token.with_suffix(".exit").exists():
        if time.monotonic() > deadline:
            raise TimeoutError("Host Slurm bridge did not respond within 120 seconds")
        time.sleep(0.2)
    code = int(token.with_suffix(".exit").read_text())
    output = token.with_suffix(".out").read_text()
    error = token.with_suffix(".err").read_text()
    for suffix in (".out", ".err", ".exit"):
        token.with_suffix(suffix).unlink()
    if code:
        raise RuntimeError(f"Host {operation} failed ({code}): {error}")
    return output


def slurm_status(output, job_id):
    for line in output.splitlines():
        fields = line.split("|")
        if len(fields) >= 3 and fields[0] == job_id:
            state = fields[1].split()[0].rstrip("+")
            if state == "COMPLETED":
                return "success" if fields[2] == "0:0" else "failed"
            if state in {
                "BOOT_FAIL",
                "CANCELLED",
                "DEADLINE",
                "FAILED",
                "NODE_FAIL",
                "OUT_OF_MEMORY",
                "PREEMPTED",
                "REVOKED",
                "TIMEOUT",
            }:
                return "failed"
    # Accounting can lag job submission and completion.
    return "running"


def main():
    operation, *arguments = sys.argv[1:]
    if operation == "submit":
        (jobscript,) = arguments
        properties = read_job_properties(jobscript)
        if properties.get("rule") == "solve_network":
            raise ValueError("Build executor must never submit solve_network")
        threads = int(properties.get("threads", 1))
        if not 1 <= threads <= 4:
            raise ValueError(f"Unexpected build thread request: {threads}")
        print(
            request("submit", [threads, Path(jobscript).resolve()])
            .strip()
            .split(";")[0]
        )
    elif operation in ("status", "cancel"):
        if not arguments or not all(job.isdecimal() for job in arguments):
            raise ValueError("Numeric Slurm job IDs required")
        output = request(operation, arguments)
        if operation == "status":
            print(slurm_status(output, arguments[0]))
    else:
        raise ValueError(f"Unknown operation {operation}")


if __name__ == "__main__":
    main()
