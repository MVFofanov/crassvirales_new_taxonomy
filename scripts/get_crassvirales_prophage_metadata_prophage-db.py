#!/usr/bin/env python3

import argparse
from pathlib import Path

import pandas as pd


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="Filter Excel table for rows where 'lineage_genomad' contains 'Crassvirales'."
    )
    parser.add_argument("--input_excel", type=Path, required=True, help="Path to the input Excel file (.xlsx)")
    parser.add_argument(
        "--output_tsv",
        type=Path,
        required=False,
        help="Path to output .tsv file (default: input name + _crassvirales.tsv)",
    )
    parser.add_argument(
        "--sheet_name", type=str, default="Table 1", help="Excel sheet name to read (default: 'Table 1')"
    )
    return parser.parse_args()


def filter_crassvirales(excel_path: Path, sheet_name: str = "Table 1") -> pd.DataFrame:
    try:
        df = pd.read_excel(excel_path, sheet_name=sheet_name)
    except Exception as e:
        raise ValueError(f"Error reading sheet '{sheet_name}' from {excel_path}: {e}")

    if "lineage_genomad" not in df.columns:
        raise KeyError("Column 'lineage_genomad' not found in the sheet.")

    filtered = df[df["lineage_genomad"].str.contains("Crassvirales", na=False)]
    return filtered


def save_dataframe_as_tsv(df: pd.DataFrame, output_path: Path) -> None:
    df.to_csv(output_path, sep="\t", index=False)


def main(input_excel: Path, output_tsv: Path | None = None, sheet_name: str = "Table 1") -> None:
    if not input_excel.exists():
        raise FileNotFoundError(f"Input file does not exist: {input_excel}")

    filtered_df = filter_crassvirales(input_excel, sheet_name)

    if output_tsv is None:
        output_tsv = input_excel.with_name(input_excel.stem + "_crassvirales.tsv")

    save_dataframe_as_tsv(filtered_df, output_tsv)
    print(f"Filtered {len(filtered_df)} rows and saved to: {output_tsv}")


if __name__ == "__main__":
    args = parse_args()
    main(args.input_excel, args.output_tsv, args.sheet_name)
