"""Genome-level annotation statistics. Run from the repo root: python -m genome_analyzer.main {stats,bgc,prokka} -h"""
import argparse
import logging
import os
import sys

from . import bgc
from . import plots
from . import prokka_update


class StreamToLogger:
    """File-like object that sends print() output to a logger."""

    def __init__(self, logger, level):
        self.logger = logger
        self.level = level

    def write(self, buf):
        for line in buf.rstrip().splitlines():
            self.logger.log(self.level, line.rstrip())

    def flush(self):
        pass


def setup_tool_logging(tool_name, output_dir):
    """Log to <output_dir>/<tool_name>.log and the console; redirect stdout/stderr to the log."""
    os.makedirs(output_dir, exist_ok=True)
    log_file = os.path.join(output_dir, f"{tool_name}.log")
    logging.basicConfig(
        level=logging.INFO,
        format='%(asctime)s - %(name)s - %(levelname)s - %(message)s',
        handlers=[logging.FileHandler(log_file, mode='w'), logging.StreamHandler(sys.stdout)],
        force=True
    )
    sys.stdout = StreamToLogger(logging.getLogger(tool_name), logging.INFO)
    sys.stderr = StreamToLogger(logging.getLogger(f"{tool_name}_ERROR"), logging.ERROR)
    logging.info(f"--- Starting {tool_name} Analysis ---")
    logging.info(f"Log file location: {log_file}")


def main():
    parser = argparse.ArgumentParser(description="Genome-level annotation statistics")
    subparsers = parser.add_subparsers(dest="command", required=True)

    p_stats = subparsers.add_parser("stats", help="Statistics and plots from the per-genome summary CSV")
    p_stats.add_argument("-i", "--input", required=True, help="Per-genome summary CSV or Excel file")
    p_stats.add_argument("-o", "--output", default="plots", help="Output directory")
    p_stats.add_argument("-g", "--gff-dir", help="Directory of GFF files, for contig counts and genome size")

    p_bgc = subparsers.add_parser("bgc", help="antiSMASH vs GECCO cluster comparison from a directory of GFF files")
    p_bgc.add_argument("-d", "--dir", required=True, help="Directory containing .gff or .gff3 files")
    p_bgc.add_argument("-o", "--output", default="bgc_results", help="Output directory")
    p_bgc.add_argument("--skip-plots", action="store_true", help="Skip statistics and plots")

    p_prokka = subparsers.add_parser("prokka", help="Add Prokka metrics to the summary CSV (edits it in place)")
    p_prokka.add_argument("--csv_file", required=True, help="Path to the summary CSV")
    p_prokka.add_argument("--prokka_gff_dir", required=True, help="Directory containing Prokka GFF files")
    p_prokka.add_argument("--prokka_gbk_dir", required=True, help="Directory containing Prokka GBK files")

    args = parser.parse_args()

    if args.command == "stats":
        setup_tool_logging("stats", args.output)
        plots.cmd_plot(args)
    elif args.command == "bgc":
        setup_tool_logging("bgc", args.output)
        bgc.cmd_plot(args)
    elif args.command == "prokka":
        setup_tool_logging("prokka", os.path.dirname(args.csv_file) or ".")
        prokka_update.update_all_genomes(args.csv_file, args.prokka_gff_dir, args.prokka_gbk_dir)


if __name__ == "__main__":
    main()
