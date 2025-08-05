#!/usr/bin/env python3

import argparse
import pandas as pd
from typing import List, Tuple

def parse_arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Cluster tblastx hits into candidate prophage regions.")
    parser.add_argument("--blast_file", required=True, help="BLAST output in outfmt 6 format.")
    parser.add_argument("--output_bed", required=True, help="Output BED-like file with merged coordinates.")
    parser.add_argument("--merge_distance", type=int, default=10000, help="Maximum distance to merge nearby hits.")
    return parser.parse_args()

def merge_intervals(intervals: List[Tuple[int, int]], max_distance: int = 10000) -> List[Tuple[int, int]]:
    intervals.sort()
    merged = []
    current_start, current_end = intervals[0]

    for start, end in intervals[1:]:
        if start <= current_end + max_distance:
            current_end = max(current_end, end)
        else:
            merged.append((current_start, current_end))
            current_start, current_end = start, end
    merged.append((current_start, current_end))
    return merged

def main():
    args = parse_arguments()
    
    # BLAST columns for outfmt 6
    cols = [
        "query", "subject", "pident", "length", "mismatch", "gapopen",
        "qstart", "qend", "sstart", "send", "evalue", "bitscore"
    ]
    
    df = pd.read_csv(args.blast_file, sep="\t", names=cols)
    df["s_min"] = df[["sstart", "send"]].min(axis=1)
    df["s_max"] = df[["sstart", "send"]].max(axis=1)
    
    records = []
    for subject, group in df.groupby("subject"):
        intervals = list(zip(group["s_min"], group["s_max"]))
        merged = merge_intervals(intervals, max_distance=args.merge_distance)
        for start, end in merged:
            records.append((subject, int(start), int(end)))
    
    out_df = pd.DataFrame(records, columns=["contig", "start", "end"])
    out_df.to_csv(args.output_bed, sep="\t", index=False, header=False)

if __name__ == "__main__":
    main()
