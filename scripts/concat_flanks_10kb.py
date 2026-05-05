#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Concatenate 10kb flanking .fna fragments (left + right) for each prophage directory
into a single FASTA file saved inside that directory as:
{prophage_name}_flanking_10kb_all_{N}.fna

Directory layout (example):
root/
  DASYBV010000020.1&provirus_269058_362578/
    DASYBV010000020.1&provirus_269058_362578_flanking_10kb_left_169058_179057.fna
    DASYBV010000020.1&provirus_269058_362578_flanking_10kb_left_179058_189057.fna
    ...
    DASYBV010000020.1&provirus_269058_362578_flanking_10kb_right_362579_372578.fna
    ...

Result:
  DASYBV010000020.1&provirus_269058_362578_flanking_10kb_all_{N}.fna
"""

from __future__ import annotations
import argparse
import logging
import re
from pathlib import Path
from typing import Iterable, List, Tuple

# Match files like:
# {prophage_name}_flanking_10kb_{left|right}_{start}_{end}.fna
FNAME_RE = re.compile(
    r"^(?P<name>.+?)_flanking_(?P<size>\d+)kb_(?P<side>left|right)_(?P<start>\d+)_(?P<end>\d+)\.fna$"
)

def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(
        description="Concatenate 10kb flanking .fna files (left+right) per prophage directory."
    )
    p.add_argument("--root", type=Path, required=True,
                   help="Root directory that contains prophage subdirectories.")
    p.add_argument("--size_kb", type=int, default=10,
                   help="Flank size in kb to include (default: 10).")
    p.add_argument("--dry_run", action="store_true",
                   help="Show what would be done without writing files.")
    p.add_argument("--overwrite", action="store_true",
                   help="Overwrite existing combined file if present.")
    p.add_argument("--verbose", "-v", action="count", default=0,
                   help="Increase log verbosity (-v, -vv).")
    return p.parse_args()

def iter_prophage_dirs(root: Path) -> Iterable[Path]:
    for p in sorted(root.iterdir()):
        if p.is_dir():
            yield p

def collect_flank_files(prophage_dir: Path, size_kb: int) -> List[Tuple[str, int, int, Path]]:
    """
    Returns list of tuples: (side, start, end, path)
    Only files matching the given flank size (kb) are included.
    """
    out: List[Tuple[str, int, int, Path]] = []
    for f in prophage_dir.glob("*.fna"):
        m = FNAME_RE.match(f.name)
        if not m:
            continue
        if int(m.group("size")) != size_kb:
            continue
        side = m.group("side")  # "left" or "right"
        start = int(m.group("start"))
        end = int(m.group("end"))
        out.append((side, start, end, f))
    # Sort deterministically: left first by start, then right by start
    # (Change order if you prefer right-first)
    out.sort(key=lambda x: (0 if x[0] == "left" else 1, x[1], x[2], x[3].name))
    return out

def read_file_text(path: Path) -> str:
    txt = path.read_text(encoding="utf-8", errors="replace")
    # Ensure trailing newline for clean concatenation
    if not txt.endswith("\n"):
        txt += "\n"
    return txt

def main() -> None:
    args = parse_args()
    log_level = logging.WARNING
    if args.verbose == 1:
        log_level = logging.INFO
    elif args.verbose >= 2:
        log_level = logging.DEBUG
    logging.basicConfig(format="%(levelname)s: %(message)s", level=log_level)

    root: Path = args.root
    if not root.exists():
        logging.error("Root directory does not exist: %s", root)
        return

    size_kb = args.size_kb

    processed = 0
    for prophage_dir in iter_prophage_dirs(root):
        prophage_name = prophage_dir.name  # e.g., DASYBV010000020.1&provirus_269058_362578
        files = collect_flank_files(prophage_dir, size_kb=size_kb)
        if not files:
            logging.info("No %dkb flanks in %s", size_kb, prophage_dir)
            continue

        # Build output file name with count
        n = len(files)
        out_name = f"{prophage_name}_flanking_{size_kb}kb_all_{n}.fna"
        out_path = prophage_dir / out_name

        if out_path.exists() and not args.overwrite:
            logging.info("Skip existing (use --overwrite to replace): %s", out_path)
            processed += 1
            continue

        # Prepare concatenated content
        parts: List[str] = []
        for side, start, end, f in files:
            logging.debug("Adding %s (%s:%d-%d)", f.name, side, start, end)
            parts.append(read_file_text(f))
        combined = "".join(parts)

        if args.dry_run:
            logging.info("[DRY RUN] Would write %s (%d files)", out_path, n)
        else:
            out_path.write_text(combined, encoding="utf-8")
            logging.info("Wrote %s (%d files)", out_path, n)

        processed += 1

    logging.info("Done. Processed %d prophage directories.", processed)

if __name__ == "__main__":
    main()
