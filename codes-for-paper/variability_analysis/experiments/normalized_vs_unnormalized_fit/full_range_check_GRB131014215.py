"""One-burst driver for the paper-facing full-range check -- GRB131014A only.
See full_range_single_burst.py for the shared logic and full_range_all_bursts_check.py's docstring for
the full rationale. Split out 2026-09-24 so each burst can run (and be killed/monitored) independently."""
from full_range_single_burst import run_one

if __name__ == "__main__":
    run_one("GRB131014215")
