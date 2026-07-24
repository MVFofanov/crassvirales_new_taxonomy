#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Concatenate flanking FASTA fragments for each prophage directory.

For every prophage, create:

1. All left and right fragments:
   {prophage_name}_flanking_{size_kb}kb_all_{N}.fna

2. Left fragments only:
   {prophage_name}_flanking_{size_kb}kb_left_{N}.fna

3. Right fragments only:
   {prophage_name}_flanking_{size_kb}kb_right_{N}.fna

By default, output files are written inside each prophage directory.
Use --outdir to collect all output files in one directory.
"""

from __future__ import annotations

import argparse
import logging
import re
from pathlib import Path
from typing import Iterable, List, Tuple


FlankFile = Tuple[str, int, int, Path]


# Match:
# {prophage_name}_flanking_10kb_{left|right}_{start}_{end}.fna
FNAME_RE = re.compile(
    r"^(?P<name>.+?)"
    r"_flanking_(?P<size>\d+)kb"
    r"_(?P<side>left|right)"
    r"_(?P<start>\d+)"
    r"_(?P<end>\d+)"
    r"\.fna$"
)


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Concatenate flanking FASTA files per prophage directory, "
            "creating combined, left-only, and right-only multi-FASTA files."
        )
    )

    parser.add_argument(
        "--root",
        type=Path,
        required=True,
        help="Root directory containing prophage subdirectories.",
    )

    parser.add_argument(
        "--outdir",
        type=Path,
        default=None,
        help=(
            "Directory for combined FASTA files. "
            "Default: write inside each prophage directory."
        ),
    )

    parser.add_argument(
        "--size_kb",
        type=int,
        default=10,
        help="Flank fragment size in kb to include (default: 10).",
    )

    parser.add_argument(
        "--dry_run",
        action="store_true",
        help="Show what would be written without creating files.",
    )

    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="Overwrite existing output files.",
    )

    parser.add_argument(
        "--verbose",
        "-v",
        action="count",
        default=0,
        help="Increase log verbosity: -v for INFO, -vv for DEBUG.",
    )

    return parser.parse_args()


def iter_prophage_dirs(root: Path) -> Iterable[Path]:
    """Yield immediate subdirectories of the root directory."""

    for path in sorted(root.iterdir()):
        if path.is_dir():
            yield path


def collect_flank_files(
    prophage_dir: Path,
    size_kb: int,
) -> List[FlankFile]:
    """
    Collect flank files matching the requested fragment size.

    Returns tuples containing:

        side, start, end, path
    """

    files: List[FlankFile] = []

    for path in prophage_dir.glob("*.fna"):
        match = FNAME_RE.match(path.name)

        if match is None:
            continue

        if int(match.group("size")) != size_kb:
            continue

        side = match.group("side")
        start = int(match.group("start"))
        end = int(match.group("end"))

        files.append((side, start, end, path))

    # Genomic coordinate order:
    # left fragments first, then right fragments.
    files.sort(
        key=lambda item: (
            0 if item[0] == "left" else 1,
            item[1],
            item[2],
            item[3].name,
        )
    )

    return files


def read_file_text(path: Path) -> str:
    """Read a FASTA file and ensure that it ends with a newline."""

    text = path.read_text(encoding="utf-8", errors="replace")

    if not text.endswith("\n"):
        text += "\n"

    return text


def write_combined_fasta(
    files: List[FlankFile],
    output_path: Path,
    overwrite: bool,
    dry_run: bool,
    label: str,
) -> bool:
    """
    Concatenate FASTA records and write them to one multi-FASTA file.

    Returns True when the file was written or would be written in dry-run mode.
    Returns False when an existing file was skipped.
    """

    if output_path.exists() and not overwrite:
        logging.info(
            "Skipping existing %s file: %s",
            label,
            output_path,
        )
        return False

    for side, start, end, path in files:
        logging.debug(
            "Adding to %s: %s (%s:%d-%d)",
            label,
            path.name,
            side,
            start,
            end,
        )

    combined = "".join(read_file_text(path) for _, _, _, path in files)

    if dry_run:
        logging.info(
            "[DRY RUN] Would write %s using %d fragments: %s",
            label,
            len(files),
            output_path,
        )
    else:
        output_path.parent.mkdir(parents=True, exist_ok=True)
        output_path.write_text(combined, encoding="utf-8")

        logging.info(
            "Wrote %s using %d fragments: %s",
            label,
            len(files),
            output_path,
        )

    return True


def main() -> None:
    args = parse_args()

    log_level = logging.WARNING

    if args.verbose == 1:
        log_level = logging.INFO
    elif args.verbose >= 2:
        log_level = logging.DEBUG

    logging.basicConfig(
        format="%(levelname)s: %(message)s",
        level=log_level,
    )

    root: Path = args.root

    if not root.is_dir():
        raise NotADirectoryError(
            f"Root directory does not exist or is not a directory: {root}"
        )

    if args.size_kb <= 0:
        raise ValueError("--size_kb must be greater than zero.")

    if args.outdir is not None and not args.dry_run:
        args.outdir.mkdir(parents=True, exist_ok=True)

    processed = 0
    outputs_written = 0

    for prophage_dir in iter_prophage_dirs(root):
        prophage_name = prophage_dir.name

        files = collect_flank_files(
            prophage_dir,
            size_kb=args.size_kb,
        )

        if not files:
            logging.info(
                "No %d kb flank files found in %s",
                args.size_kb,
                prophage_dir,
            )
            continue

        left_files = [
            item for item in files
            if item[0] == "left"
        ]

        right_files = [
            item for item in files
            if item[0] == "right"
        ]

        output_directory = (
            args.outdir
            if args.outdir is not None
            else prophage_dir
        )

        output_jobs = [
            (
                "all flanks",
                files,
                output_directory
                / (
                    f"{prophage_name}_flanking_"
                    f"{args.size_kb}kb_all_{len(files)}.fna"
                ),
            )
        ]

        if left_files:
            output_jobs.append(
                (
                    "left flanks",
                    left_files,
                    output_directory
                    / (
                        f"{prophage_name}_flanking_"
                        f"{args.size_kb}kb_left_{len(left_files)}.fna"
                    ),
                )
            )
        else:
            logging.warning(
                "No left flank fragments found for %s",
                prophage_name,
            )

        if right_files:
            output_jobs.append(
                (
                    "right flanks",
                    right_files,
                    output_directory
                    / (
                        f"{prophage_name}_flanking_"
                        f"{args.size_kb}kb_right_{len(right_files)}.fna"
                    ),
                )
            )
        else:
            logging.warning(
                "No right flank fragments found for %s",
                prophage_name,
            )

        for label, selected_files, output_path in output_jobs:
            written = write_combined_fasta(
                files=selected_files,
                output_path=output_path,
                overwrite=args.overwrite,
                dry_run=args.dry_run,
                label=label,
            )

            if written:
                outputs_written += 1

        logging.info(
            "Processed %s: %d left, %d right, %d total fragments.",
            prophage_name,
            len(left_files),
            len(right_files),
            len(files),
        )

        processed += 1

    logging.warning(
        "Done. Processed %d prophage directories; created or planned %d output files.",
        processed,
        outputs_written,
    )


if __name__ == "__main__":
    main()