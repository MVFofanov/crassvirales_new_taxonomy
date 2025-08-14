#!/usr/bin/env python3

from pathlib import Path

import pandas as pd


def read_length_table(path: Path) -> pd.DataFrame:
    """Reads a 2-column table with contig_id and length."""
    return pd.read_csv(path, sep="\t", names=["id", "length"], header=0, dtype={"id": str})


def read_prophage_lengths(path: Path) -> pd.DataFrame:
    """Reads a 2-column table with prophage_id and length."""
    return pd.read_csv(path, sep="\t", names=["prophage_id", "prophage_length"], header=0, dtype={"prophage_id": str})


def read_prophage_table(path: Path) -> pd.DataFrame:
    """Reads the main prophage coordinates table."""
    return pd.read_csv(
        path,
        # sep="\t",
        delim_whitespace=True,  # ✅ auto-detects spaces or tabs
        names=["prophage_id", "contig_id", "start", "end", "predicted_length"],
        header=None,
        dtype={"contig_id": str},
    )


def add_lengths(prophage_df: pd.DataFrame, contig_lengths_df: pd.DataFrame) -> pd.DataFrame:
    """Adds contig lengths to the prophage dataframe."""
    merged = prophage_df.merge(
        contig_lengths_df.rename(columns={"id": "contig_id", "length": "contig_length"}), on="contig_id", how="left"
    )
    return merged


def save_table(df: pd.DataFrame, output_path: Path) -> None:
    """Saves the DataFrame to a TSV file."""
    df.to_csv(output_path, sep="\t", index=False)


def main():
    wd = "/mnt/c/crassvirales/crassvirales_new_taxonomy/crassvirales_prophages"
    blast_dir = f"{wd}/blast_prophages_vs_ncbi_and_gtdb"
    prophage_table_path = Path(f"{blast_dir}/prophage_alignments_ncbi_and_gtdb_samtools_coordinates_edited.txt")
    contig_lengths_path = Path(f"{blast_dir}/TerL_crassvirales_genome_ids_ncbi_and_gtdb_genome_size.tsv")
    prophage_lengths_path = Path(f"{wd}/crassvirales_prophages_ncbi_and_prophage-db_genome_size.tsv")
    output_path = Path(f"{blast_dir}/prophage_alignments_ncbi_and_gtdb_samtools_coordinates_edited_with_lengths.tsv")

    # Read input tables
    prophage_df = read_prophage_table(prophage_table_path)
    contig_lengths_df = read_length_table(contig_lengths_path)
    prophage_lengths_df = read_prophage_lengths(prophage_lengths_path)

    # Merge contig lengths
    result_df = add_lengths(prophage_df, contig_lengths_df)

    # Merge prophage lengths
    result_df = result_df.merge(prophage_lengths_df, on="prophage_id", how="left")

    # Reorder columns
    columns = ["prophage_id", "contig_id", "start", "end", "predicted_length", "prophage_length", "contig_length"]
    result_df = result_df.reindex(columns=columns)

    # Save output
    save_table(result_df, output_path)
    print(f"Saved merged table with lengths to: {output_path}")


if __name__ == "__main__":
    main()
