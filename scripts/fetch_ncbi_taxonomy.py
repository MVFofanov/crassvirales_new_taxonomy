#!/usr/bin/env python3

import argparse
import time
from pathlib import Path

import pandas as pd
from Bio import Entrez

Entrez.email = "mikhail.v.fofanov@gmail.com"  # Change this to your email


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Fetch taxonomy lineages from NCBI for a list of IDs.")
    parser.add_argument(
        "--id_type", choices=["assembly", "genome", "nucleotide"], required=True, help="Type of NCBI ID"
    )
    parser.add_argument("--id_list", type=Path, required=True, help="Path to file with one ID per line")
    parser.add_argument("--taxonomy_output", type=Path, required=True, help="Output TSV file with taxonomy info")
    return parser.parse_args()


def read_id_list(path: Path) -> list[str]:
    with path.open() as f:
        return [line.strip() for line in f if line.strip()]


def fetch_taxonomy_for_id(id_type: str, accession: str) -> tuple[str, str]:
    try:
        if id_type == "assembly":
            search_handle = Entrez.esearch(db="assembly", term=accession)
            search_results = Entrez.read(search_handle)
            search_handle.close()

            # print(search_results["IdList"][0])
            # print(search_results["IdList"])
            # print(search_results)

            if not search_results["IdList"]:
                return accession, "NA"

            uid = search_results["IdList"][0]

            summary_handle = Entrez.esummary(db="assembly", id=uid, report="full")
            summary = Entrez.read(summary_handle, validate=False)
            summary_handle.close()

            taxid = summary["DocumentSummarySet"]["DocumentSummary"][0]["Taxid"]
            tax_handle = Entrez.efetch(db="taxonomy", id=taxid)
            tax_data = Entrez.read(tax_handle)
            tax_handle.close()

            lineage = tax_data[0]["Lineage"]
            print(f"{lineage=}")
            # print(summary['DocumentSummarySet'])
            # print(summary['DocumentSummarySet']['DocumentSummary'])
            # print(summary['DocumentSummarySet']['DocumentSummary'][0])
            return accession, lineage

        elif id_type == "genome":
            search_handle = Entrez.esearch(db="nuccore", term=f"{accession}[ACCN]")
            search_results = Entrez.read(search_handle)
            search_handle.close()

            if not search_results["IdList"]:
                return accession, "NA"

            uid = search_results["IdList"][0]
            summary_handle = Entrez.esummary(db="nuccore", id=uid)
            summary = Entrez.read(summary_handle, validate=False)
            summary_handle.close()

            taxid = summary[0]["TaxId"]
            tax_handle = Entrez.efetch(db="taxonomy", id=taxid)
            tax_data = Entrez.read(tax_handle)
            tax_handle.close()

            lineage = tax_data[0]["Lineage"]
            return accession, lineage

        elif id_type == "nucleotide":
            handle = Entrez.efetch(db="nucleotide", id=accession, rettype="gb", retmode="xml")
            records = Entrez.read(handle)
            handle.close()

            taxid = records[0]["GBSeq_taxonomy"]
            return accession, taxid

    except Exception as e:
        print(f"[ERROR] Failed to fetch taxonomy for {accession}: {e}")
        return accession, "NA"


def fetch_taxonomies(id_type: str, id_list: list[str]) -> list[tuple[str, str]]:
    results = []
    for acc in id_list:
        acc_clean = acc.strip()
        print(f"[INFO] Fetching taxonomy for {acc_clean}...")
        result = fetch_taxonomy_for_id(id_type, acc_clean)
        results.append(result)
        time.sleep(0.34)  # Respect NCBI rate limit (3 requests/sec)
    return results


def save_taxonomy_table(results: list[tuple[str, str]], output_path: Path) -> None:
    df = pd.DataFrame(results, columns=["accession", "taxonomy"])
    df.to_csv(output_path, sep="\t", index=False)


def main():
    args = parse_args()
    id_list = read_id_list(args.id_list)
    taxonomy_data = fetch_taxonomies(args.id_type, id_list)
    save_taxonomy_table(taxonomy_data, args.taxonomy_output)
    print(f"[DONE] Taxonomy table saved to {args.taxonomy_output}")


if __name__ == "__main__":
    main()
