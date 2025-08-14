#!/usr/bin/env python3

import argparse
from pathlib import Path

from Bio import SeqIO


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Extract sequences by ID from a multi-FASTA file.")
    parser.add_argument("--input_fasta", required=True, type=Path, help="Path to multi-FASTA input file.")
    parser.add_argument(
        "--id_file", required=True, type=Path, help="Path to text file with sequence IDs (one per line)."
    )
    parser.add_argument("--output_fasta", required=True, type=Path, help="Path to output FASTA file.")
    return parser.parse_args()


def load_ids(id_file: Path) -> set[str]:
    with id_file.open() as f:
        return set(line.strip() for line in f if line.strip())


def extract_sequences(input_fasta: Path, target_ids: set[str]) -> list:
    return [
        record
        for record in SeqIO.parse(str(input_fasta), "fasta")
        if record.id in target_ids or record.description in target_ids
    ]


def write_fasta(sequences: list, output_fasta: Path) -> None:
    if sequences:
        SeqIO.write(sequences, str(output_fasta), "fasta")
        print(f"✅ Wrote {len(sequences)} sequences to {output_fasta}")
    else:
        print("⚠️ No matching sequences found.")


def main() -> None:
    args = parse_args()
    ids = load_ids(args.id_file)
    sequences = extract_sequences(args.input_fasta, ids)
    write_fasta(sequences, args.output_fasta)


if __name__ == "__main__":
    main()
