"""Copy every figure and generated table this paper uses from GRBResearchWork into
GRBResearchPaper, per `figure_registry.yaml` and `table_registry.yaml`.

Root-level, same reasoning as `seed_table_to_latex.py`: this spans every topic folder's
own output rather than belonging to one of them, and it's paper-output plumbing, not
analysis code.

This only copies files that already exist -- it does not regenerate anything. Run each
topic's own script(s) first if a source is stale or missing; this script will say exactly
which file is missing rather than silently skipping it.

Usage: PYTHONPATH=../src .venv/bin/python codes-for-paper/sync_paper_assets.py
"""

import shutil
import sys
from pathlib import Path

import yaml

HERE = Path(__file__).parent
WORK_ROOT = HERE.parent
PAPER_ROOT = WORK_ROOT.parent / "GRBResearchPaper"

REGISTRIES = [
    ("figure", HERE / "figure_registry.yaml"),
    ("table", HERE / "table_registry.yaml"),
]


def sync_one(kind: str, registry_path: Path) -> tuple[int, int]:
    with open(registry_path) as f:
        registry = yaml.safe_load(f)

    copied = 0
    missing = []
    for slug, entry in registry.get("active", {}).items():
        src = WORK_ROOT / entry["source"]
        dest = PAPER_ROOT / entry["dest"]
        if not src.is_file():
            missing.append((slug, src))
            continue
        dest.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(src, dest)
        print(f"  [{kind}] {slug}: {entry['source']} -> {entry['dest']}")
        copied += 1

    for slug, src in missing:
        print(f"  [{kind}] MISSING, not copied: {slug} ({src})", file=sys.stderr)

    n_unused = len(registry.get("unused", {}))
    n_no_source = len(registry.get("no_known_source", {}))
    if n_unused or n_no_source:
        print(f"  [{kind}] {n_unused} unused, {n_no_source} no_known_source entry/entries skipped by design")

    return copied, len(missing)


def main():
    if not PAPER_ROOT.is_dir():
        raise SystemExit(f"GRBResearchPaper not found as a sibling of GRBResearchWork: {PAPER_ROOT}")

    total_copied = 0
    total_missing = 0
    for kind, registry_path in REGISTRIES:
        print(f"=== {kind}s ({registry_path.name}) ===")
        copied, missing = sync_one(kind, registry_path)
        total_copied += copied
        total_missing += missing

    print(f"\n{total_copied} file(s) copied into {PAPER_ROOT}.")
    if total_missing:
        print(f"{total_missing} active entry/entries had a missing source file -- see stderr above.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
