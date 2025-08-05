#!/usr/bin/env python3

import argparse
from pathlib import Path


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Remove last two fields from pipe-delimited TerL labels.")
    parser.add_argument("--input_file", required=True, type=Path, help="Input text file with TerL labels.")
    parser.add_argument("--output_file", required=True, type=Path, help="Output file with genome IDs (without last two fields).")
    return parser.parse_args()


def strip_last_two_fields(line: str, delimiter: str = "|") -> str:
    parts = line.strip().split(delimiter)
    if len(parts) <= 2:
        return ""  # Or raise warning/error if needed
    return delimiter.join(parts[:-2])


def process_file(input_file: Path, output_file: Path) -> None:
    with input_file.open("r") as fin, output_file.open("w") as fout:
        count = 0
        for line in fin:
            result = strip_last_two_fields(line)
            if result:
                fout.write(result + "\n")
                count += 1
        print(f"✅ Processed {count} lines. Output saved to: {output_file}")


def main() -> None:
    args = parse_args()
    process_file(args.input_file, args.output_file)


if __name__ == "__main__":
    main()
