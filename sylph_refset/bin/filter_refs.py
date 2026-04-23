#!/usr/bin/env python3

import argparse
import csv
import re
from pathlib import Path


def parse_args():
    parser = argparse.ArgumentParser(
        description="Filter a combined sylph TSV report using ANI and effective coverage thresholds."
    )
    parser.add_argument("--input", type=Path, required=True, help="Combined sylph report TSV.")
    parser.add_argument(
        "--taxonomy-data",
        type=Path,
        help="Two-column taxonomy TSV mapping genome IDs to GTDB taxonomy strings.",
    )
    parser.add_argument(
        "--sylph-method",
        choices=["query", "profile"],
        required=True,
        help="Sylph mode used to produce the report.",
    )
    parser.add_argument("--ani", type=float, required=True, help="Minimum Adjusted_ANI threshold.")
    parser.add_argument("--cov", type=float, required=True, help="Minimum Eff_cov threshold.")
    parser.add_argument("--out-report", type=Path, required=True, help="Filtered TSV output path.")
    parser.add_argument("--out-summary", type=Path, help="Optional summary TSV output path.")
    parser.add_argument(
        "--out-references",
        type=Path,
        help="Optional file containing unique Genome_file entries from the filtered report.",
    )
    return parser.parse_args()


def parse_float(value, column, row_number):
    try:
        return float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"Row {row_number} has a non-numeric value in '{column}': {value!r}") from exc


def extract_genome_id(genome_file):
    match = re.search(r"/([^/]+)_genomic\.fna\.gz$", genome_file or "")
    if match:
        return match.group(1)
    return Path(genome_file).name if genome_file else ""


def resolve_column(fieldnames, candidates, label):
    for candidate in candidates:
        if candidate in fieldnames:
            return candidate
    raise ValueError(
        f"Input report is missing a supported {label} column. Tried: {', '.join(candidates)}"
    )


def get_filter_columns(method, fieldnames):
    if method == "query":
        ani_column = resolve_column(fieldnames, ["Adjusted_ANI"], "ANI")
        cov_column = resolve_column(fieldnames, ["Eff_cov"], "coverage")
        return ani_column, cov_column

    ani_column = resolve_column(fieldnames, ["Naive_ANI", "Naive ANI"], "ANI")
    cov_column = resolve_column(fieldnames, ["Eff_cov", "eff cov", "Eff cov", "eff_cov"], "coverage")
    return ani_column, cov_column


def parse_species_name(taxonomy_string):
    match = re.search(r"s__([^|;]*)", taxonomy_string or "")
    if match and match.group(1):
        return match.group(1)
    return "unknown_species"


def load_species_by_genome_id(path):
    species_by_genome_id = {}
    if not path:
        return species_by_genome_id

    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        for row in reader:
            if len(row) < 2:
                continue
            species_by_genome_id[row[0]] = parse_species_name(row[1])
    return species_by_genome_id


def main():
    args = parse_args()
    species_by_genome_id = load_species_by_genome_id(args.taxonomy_data)

    with args.input.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        fieldnames = reader.fieldnames or []
        ani_column, cov_column = get_filter_columns(args.sylph_method, fieldnames)
        required_columns = {ani_column, cov_column, "Genome_file"}
        missing_columns = sorted(required_columns.difference(fieldnames))
        if missing_columns:
            raise ValueError(
                f"Input report is missing required columns: {', '.join(missing_columns)}"
            )

        filtered_rows = []
        unique_references = []
        seen_references = set()
        all_genome_ids = set()
        retained_genome_ids = set()
        total_rows = 0

        for row_number, row in enumerate(reader, start=2):
            total_rows += 1
            genome_id = extract_genome_id(row["Genome_file"])
            if genome_id:
                all_genome_ids.add(genome_id)
            ani = parse_float(row.get(ani_column), ani_column, row_number)
            cov = parse_float(row.get(cov_column), cov_column, row_number)
            if ani < args.ani or cov < args.cov:
                continue

            filtered_rows.append(row)
            reference = row["Genome_file"]
            if reference and reference not in seen_references:
                seen_references.add(reference)
                unique_references.append(reference)
            if genome_id:
                retained_genome_ids.add(genome_id)

    args.out_report.parent.mkdir(parents=True, exist_ok=True)
    with args.out_report.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t")
        writer.writeheader()
        writer.writerows(filtered_rows)

    if args.out_summary:
        removed_genome_ids = sorted(all_genome_ids.difference(retained_genome_ids))
        removed_by_species = {}
        for genome_id in removed_genome_ids:
            species_name = species_by_genome_id.get(genome_id, "unknown_species")
            removed_by_species.setdefault(species_name, []).append(genome_id)

        args.out_summary.parent.mkdir(parents=True, exist_ok=True)
        with args.out_summary.open("w", newline="") as handle:
            writer = csv.DictWriter(
                handle,
                fieldnames=["species_name", "removed_reference_count", "removed_genome_ids"],
                delimiter="\t",
            )
            writer.writeheader()
            for species_name in sorted(removed_by_species):
                genome_ids = sorted(set(removed_by_species[species_name]))
                writer.writerow(
                    {
                        "species_name": species_name,
                        "removed_reference_count": len(genome_ids),
                        "removed_genome_ids": ",".join(genome_ids),
                    }
                )

    if args.out_references:
        args.out_references.parent.mkdir(parents=True, exist_ok=True)
        with args.out_references.open("w", newline="") as handle:
            for reference in unique_references:
                handle.write(f"{reference}\n")


if __name__ == "__main__":
    main()
