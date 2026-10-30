"""Independent wall-time supervision, including Julia startup and final auditing."""
import json
import os
from pathlib import Path
import signal
import subprocess
import sys
import time


def supervise(command, out, budget):
    out = Path(out)
    log = out.with_name(out.name + ".log")
    if out.exists() or log.exists():
        raise ValueError("output or log already exists")
    if not 0 < budget < float("inf"):
        raise ValueError("budget must be finite and positive")
    out.parent.mkdir(parents=True, exist_ok=True)
    env = dict(os.environ, AD_DEADLINE=str(time.time() + budget))
    start = time.monotonic()
    with log.open("x") as stream:
        child = subprocess.Popen(command, stdout=stream, stderr=subprocess.STDOUT,
                                 env=env, start_new_session=True)
        timed_out = False
        try:
            child.wait(timeout=budget)
        except subprocess.TimeoutExpired:
            timed_out = True
            os.killpg(child.pid, signal.SIGTERM)
            try:
                child.wait(timeout=5)
            except subprocess.TimeoutExpired:
                os.killpg(child.pid, signal.SIGKILL)
                child.wait()
    out.mkdir(exist_ok=True)
    checkpoints = sorted(out.glob("step[0-9][0-9][0-9][0-9].json"))
    latest = checkpoints[-1].name if checkpoints else "initial.json" if (out / "initial.json").exists() else None
    supervision = dict(command=command, budget=budget, seconds=time.monotonic() - start,
                       exit=child.returncode, timed_out=timed_out, last_checkpoint=latest)
    (out / "supervision.json").write_text(json.dumps(supervision, indent=2) + "\n")
    summary_path = out / "summary.json"
    if timed_out or child.returncode != 0 or not summary_path.exists():
        result = dict(status="budget" if timed_out else "worker_failed", accepted=False,
                      **supervision)
        temporary = out / "summary-supervisor.tmp"
        temporary.write_text(json.dumps(result, indent=2) + "\n")
        temporary.replace(summary_path)
        return result
    return json.loads(summary_path.read_text())


if __name__ == "__main__":
    options = dict(arg[2:].split("=", 1) for arg in sys.argv[1:])
    root = Path(__file__).resolve().parents[2]
    entry = "campaign.jl" if options.get("campaign") == "true" else "run.jl"
    command = ["julia", "--startup-file=no", "-O1", f"--heap-size-hint={os.environ.get('AD_HEAP_HINT', '2G')}", f"--project={root / 'julia-env'}",
               str(Path(__file__).with_name(entry)), *sys.argv[1:]]
    print(json.dumps(supervise(command, options["out"], float(options.get("budget", 1800))), indent=2))
