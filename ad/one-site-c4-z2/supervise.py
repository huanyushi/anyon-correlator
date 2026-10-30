"""Run the analytic or random Julia entry with an external time budget."""
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
    progress_path = out / "progress.json"
    progress = json.loads(progress_path.read_text()) if progress_path.exists() else {}
    active = Path(progress.get("active_directory", out))
    if progress and active == out:
        attempts = sorted(path for path in out.glob("attempt*-seed*") if path.is_dir())
        if attempts:
            active = attempts[-1]
    if not active.is_dir():
        active = out
    checkpoints = sorted(active.glob("step[0-9][0-9][0-9][0-9].json"))
    latest_path = checkpoints[-1] if checkpoints else active / "initial.json"
    latest = str(latest_path.relative_to(out)) if latest_path.exists() else None
    supervision = dict(command=command, budget=budget, seconds=time.monotonic() - start,
                       exit=child.returncode, timed_out=timed_out, last_checkpoint=latest)
    (out / "supervision.json").write_text(json.dumps(supervision, indent=2) + "\n")
    summary_path = out / "summary.json"
    if timed_out or child.returncode != 0 or not summary_path.exists():
        result = dict(status="budget" if timed_out else "worker_failed", accepted=False,
                      **{key: progress[key] for key in ("attempts", "best", "active_directory") if key in progress},
                      **supervision)
        temporary = out / "summary-supervisor.tmp"
        temporary.write_text(json.dumps(result, indent=2) + "\n")
        temporary.replace(summary_path)
        return result
    return json.loads(summary_path.read_text())


if __name__ == "__main__":
    options = dict(arg[2:].split("=", 1) for arg in sys.argv[1:])
    root = Path(__file__).resolve().parents[2]
    if "resume" in options or "campaign" in options:
        raise ValueError("use --init=exact or --init=random; manual continuation is no longer an entry")
    entry = "random.jl" if options.get("init", "exact") == "random" else "run.jl"
    command = ["julia", "--startup-file=no", "-O1", f"--heap-size-hint={os.environ.get('AD_HEAP_HINT', '2G')}", f"--project={root / 'julia-env'}",
               str(Path(__file__).with_name(entry)), *sys.argv[1:]]
    print(json.dumps(supervise(command, options["out"], float(options.get("budget", 1800))), indent=2))
