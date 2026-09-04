"""Run every real pipeline script in this project, in order, as real subprocesses.

Root-level, same reasoning as `seed_table_to_latex.py`/`sync_paper_assets.py`: this spans every
topic folder rather than belonging to one of them, and it's project plumbing, not analysis code.

Built for exactly the situation that motivated it: a shared constant (N_GRID) changed, and
"which scripts need rerunning" had to be re-derived from memory. Reads `runner_registry.yaml`
for the script list (see that file for what's deliberately excluded, and why) and runs each one
exactly as a human would -- a real subprocess, `cwd` set to the script's own directory,
`PYTHONPATH` set to `src` -- since several of these scripts execute their logic at bare module
level with no `main()`/`if __name__` guard, so importing them in-process isn't a safe option.

This only *runs* scripts that already exist; it does not regenerate `runner_registry.yaml`
itself or invent new scripts to run.

Usage: .venv/bin/python codes-for-paper/runner_all.py [--only SUBSTRING] [--dry-run]
"""

import argparse
import subprocess
import sys
import time
from pathlib import Path

import yaml

HERE = Path(__file__).parent
WORK_ROOT = HERE.parent
SRC = WORK_ROOT / "src"
PYTHON = WORK_ROOT / ".venv" / "bin" / "python"
REGISTRY_PATH = HERE / "runner_registry.yaml"
LOG_DIR = HERE / "runner_all_logs"


def load_scripts():
    with open(REGISTRY_PATH) as f:
        registry = yaml.safe_load(f)
    return [entry["path"] for entry in registry["scripts"]]


def run_one(rel_path: str) -> dict:
    script_path = HERE / rel_path
    if not script_path.is_file():
        return {"path": rel_path, "status": "MISSING", "elapsed": 0.0, "log": None}

    log_path = LOG_DIR / (rel_path.replace("/", "__") + ".log")
    log_path.parent.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    with open(log_path, "w") as log_file:
        result = subprocess.run(
            [str(PYTHON), script_path.name],
            cwd=script_path.parent,
            env={"PYTHONPATH": str(SRC), "PATH": "/usr/bin:/bin"},
            stdout=log_file,
            stderr=subprocess.STDOUT,
        )
    elapsed = time.time() - t0
    status = "OK" if result.returncode == 0 else f"FAILED (exit {result.returncode})"
    return {"path": rel_path, "status": status, "elapsed": elapsed, "log": str(log_path)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--only", default=None, help="Only run scripts whose path contains this substring.")
    parser.add_argument("--dry-run", action="store_true", help="Print the script list and exit.")
    args = parser.parse_args()

    scripts = load_scripts()
    if args.only:
        scripts = [s for s in scripts if args.only in s]

    if args.dry_run:
        for s in scripts:
            print(s)
        return

    print(f"Running {len(scripts)} script(s); logs in {LOG_DIR}/\n")

    results = []
    for rel_path in scripts:
        print(f"  -> {rel_path} ...", end=" ", flush=True)
        result = run_one(rel_path)
        results.append(result)
        print(f"{result['status']}  ({result['elapsed']:.1f}s)")

    print("\n=== Summary ===")
    total = sum(r["elapsed"] for r in results)
    failed = [r for r in results if r["status"] != "OK"]
    for r in results:
        print(f"  {r['status']:<18} {r['elapsed']:>7.1f}s  {r['path']}")
    print(f"\nTotal wall time: {total:.1f}s ({total / 60:.1f} min)")

    if failed:
        print(f"\n{len(failed)} script(s) failed -- see their log files:")
        for r in failed:
            print(f"  {r['path']}: {r['log']}")
        sys.exit(1)
    else:
        print("\nAll scripts succeeded.")


if __name__ == "__main__":
    main()
