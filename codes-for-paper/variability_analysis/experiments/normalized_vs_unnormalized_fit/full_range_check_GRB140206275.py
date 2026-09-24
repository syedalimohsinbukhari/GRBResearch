"""One-burst driver for the paper-facing full-range check -- GRB140206B only.
See full_range_single_burst.py for the shared logic and full_range_all_bursts_check.py's docstring for
the full rationale. Split out 2026-09-24 so each burst can run (and be killed/monitored) independently.
This is the burst most likely to be the slow one (widest window, 7 pulses) -- run it first/alone if the
other three finish quickly."""
from full_range_single_burst import run_one

if __name__ == "__main__":
    run_one("GRB140206275")
