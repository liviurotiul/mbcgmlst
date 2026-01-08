import argparse
import sys
from pathlib import Path

from . import __version__
from .pipeline import run_pipeline

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="cgMLST allele calling using minimap2 and Ridom-style output."
    )
    parser.add_argument(
        "--alleles",
        required=True,
        help=(
            "Allele FASTA file or directory of FASTA files (one per locus). "
            "Directory inputs may contain numeric-only headers with the locus in the filename."
        ),
    )
    parser.add_argument(
        "--genome-dir",
        help="Directory of genome FASTA files (fa/fna/fasta).",
    )
    parser.add_argument(
        "--output",
        required=True,
        help="Output CSV path (Ridom-style, semicolon-delimited).",
    )
    parser.add_argument(
        "genomes",
        nargs="*",
        help="Genome FASTA files (fa/fna/fasta).",
    )
    parser.add_argument(
        "--min-identity",
        type=float,
        default=0.95,
        help="Minimum alignment identity (0-1) for allele hits (default: 0.95).",
    )
    parser.add_argument(
        "--strict",
        action="store_true",
        help="Require full-length alignments with no indels or clipping.",
    )
    parser.add_argument(
        "--allow-locus-only",
        action="store_true",
        help=(
            "Allow allele FASTA headers that contain only the locus name; "
            "allele numbers are assigned in input order."
        ),
    )
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Threads for minimap2 (default: 1).",
    )
    parser.add_argument(
        "--missing-value",
        default="NA",
        help="Value to use when no allele is called (default: NA).",
    )
    parser.add_argument(
        "--novel-fasta",
        help="Path to write novel allele sequences (default: <output>_novel_alleles.fasta).",
    )
    parser.add_argument(
        "--tmp-dir",
        action="store_true",
        help="Create a temporary working directory under the current directory (random name).",
    )
    parser.add_argument(
        "--keep-temp",
        action="store_true",
        help="Keep temporary minimap2 files instead of deleting them.",
    )
    parser.add_argument(
        "--test",
        action="store_true",
        help="Run on only one genome (the first in the provided list).",
    )
    parser.add_argument("--version", action="version", version=f"mbcgmlst {__version__}")
    return parser


def main(argv=None) -> int:
    if argv is None:
        argv = sys.argv[1:]
    if not argv:
        try:
            from .gui import main as gui_main
        except Exception as exc:
            print(f"Error: {exc}", file=sys.stderr)
            return 1
        return gui_main()

    parser = build_parser()
    args = parser.parse_args(argv)

    if not 0 < args.min_identity <= 1:
        parser.error("--min-identity must be between 0 and 1.")
    if args.threads < 1:
        parser.error("--threads must be >= 1.")

    try:
        tmp_dir = str(Path.cwd()) if args.tmp_dir else None
        progress_state = {"last_len": 0, "bar_width": 28}

        def progress_callback(current: int, total: int, sample_id: str) -> None:
            filled = int(progress_state["bar_width"] * current / total)
            bar = "#" * filled + "-" * (progress_state["bar_width"] - filled)
            message = f"[{bar}] {current}/{total} {sample_id}"
            pad = max(0, progress_state["last_len"] - len(message))
            progress_state["last_len"] = len(message)
            end = "\n" if current == total else "\r"
            print(message + (" " * pad), end=end, file=sys.stderr, flush=True)

        run_pipeline(
            alleles_path=args.alleles,
            genome_paths=args.genomes,
            genome_dir=args.genome_dir,
            output_path=args.output,
            min_identity=args.min_identity,
            strict=args.strict,
            threads=args.threads,
            missing_value=args.missing_value,
            novel_fasta=args.novel_fasta,
            tmp_dir=tmp_dir,
            keep_temp=args.keep_temp,
            allow_locus_only=args.allow_locus_only,
            max_genomes=1 if args.test else None,
            progress_callback=progress_callback,
        )
    except Exception as exc:
        print(f"Error: {exc}", file=sys.stderr)
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
