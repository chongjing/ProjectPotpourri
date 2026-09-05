#!/usr/bin/env python3
"""Merge xml/run_receipts.tsv into accessions.tsv / accessions_registered.tsv
and write data_availability_statement.md."""
from __future__ import annotations

import csv
from pathlib import Path

ROOT = Path(__file__).resolve().parent
STUDY = "PRJEB125245"
ERP = "ERP204980"
HOLD = "2027-09-05"


def main() -> None:
    recs = {}
    with (ROOT / "xml" / "run_receipts.tsv").open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            recs[row["library_name"]] = row

    acc_path = ROOT / "accessions.tsv"
    rows = list(csv.DictReader(acc_path.open(), delimiter="\t"))
    fields = rows[0].keys() if rows else []
    for row in rows:
        lib = row["library_name"]
        if lib in recs and recs[lib].get("success") == "true":
            row["run_accession"] = recs[lib]["ERR"]
            row["experiment_accession"] = recs[lib]["ERX"]
    with acc_path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=list(rows[0].keys()), delimiter="\t")
        w.writeheader()
        w.writerows(rows)

    meta = {r["library_name"]: r for r in csv.DictReader((ROOT / "metadata.tsv").open(), delimiter="\t")}
    das = ROOT / "data_availability_statement.md"
    lines = [
        "# Data Availability Statement",
        "",
        "All raw RNA-seq reads generated in this study have been deposited in the",
        "European Nucleotide Archive (ENA) at EMBL-EBI under project accession",
        f"**{STUDY}** (https://www.ebi.ac.uk/ena/browser/view/{STUDY}).",
        "Data are held private until "
        f"**{HOLD}** and can be released earlier from the Webin portal when the paper is accepted.",
        "",
        "| Library | Genotype | Pi | Inoculation | Replicate | Sample (ERS) | Run (ERR) | Experiment (ERX) |",
        "|---|---|---|---|---|---|---|---|",
    ]
    for row in rows:
        if row["library_name"] == "STUDY":
            continue
        m = meta[row["library_name"]]
        lines.append(
            f"| {row['library_name']} | {m['genotype']} | {m['phosphate']} {m['phosphate_uM']} uM | "
            f"{m['inoculation']} | {m['replicate']} | {row['sample_accession']} | "
            f"{row.get('run_accession','')} | {row.get('experiment_accession','')} |"
        )
    lines += [
        "",
        "**Accession summary:**",
        f"- BioProject: {STUDY} (study {ERP})",
        "- 32 samples, 64 paired-end FASTQ files, ~146 GB",
        "- Instrument: Illumina NovaSeq 6000, 150 bp paired-end mRNA (oligo-dT), Novogene Cambridge",
        "",
    ]
    das.write_text("\n".join(lines))
    print(f"Updated {acc_path}")
    print(f"Wrote {das}")
    ok = sum(1 for r in recs.values() if r.get("success") == "true")
    print(f"runs ok={ok}/{len(recs)}")


if __name__ == "__main__":
    main()
