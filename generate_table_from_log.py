"""
Script to generate LaTeX tables from GRB spectral analysis log files.

Generates one LaTeX table per GRB found in the log file (LogParser.generate_multiple_latex_tables,
always multi-GRB mode -- this script previously accepted a single output file and GRB name, but
neither was ever actually used: parse_log_and_generate_table defaults to multi_grb=True regardless,
which ignores both in favor of one file per burst, named from short_to_long).

Usage:
    python generate_table_from_log.py <log_file> [output_dir]

Examples:
    python generate_table_from_log.py cstat_run_20260204_172752.log
    python generate_table_from_log.py cstat_run_20260204_172752.log appendices/
"""

import sys

from src.grb_research.log_to_latex_parser import parse_log_and_generate_table


def main():
    if len(sys.argv) < 2:
        print(__doc__)
        sys.exit(1)

    log_file = sys.argv[1]
    output_dir = sys.argv[2] if len(sys.argv) > 2 else None

    try:
        parse_log_and_generate_table(log_file, output_dir)
    except Exception as e:
        print(f"Error: {e}")
        import traceback

        traceback.print_exc()
        sys.exit(1)


if __name__ == "__main__":
    main()
