#!/usr/bin/env python3
"""After Webin portal registration, put PRJEB / ERS accessions in accessions.tsv
then re-run generate_ena_pack.py so manifests switch from aliases to accessions.

accessions.tsv columns:
  library_name, sample_alias, sample_accession, study_accession, ...

The STUDY row holds the project accession. Each r2310_* row holds ERS.
"""
from __future__ import annotations

import csv
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
PATH = ROOT / "accessions.tsv"


def main() -> int:
    if len(sys.argv) < 3:
        print(
            "Usage:\n"
            "  python3 fill_accessions.py STUDY PRJEB123456\n"
            "  python3 fill_accessions.py r2310_1 ERS12345678\n"
            "  python3 fill_accessions.py --from-tsv path.tsv   # two columns: alias_or_lib  accession"
        )
        return 1
    if sys.argv[1] == "--from-tsv":
        src = Path(sys.argv[2])
        updates = {}
        study = None
        with src.open() as fh:
            for line in fh:
                if not line.strip() or line.startswith("#"):
                    continue
                a, b = line.split()[:2]
                if a.upper() in {"STUDY", "PRJEB"} or a.startswith("PRJEB"):
                    study = b if b.startswith("PRJ") else a
                else:
                    updates[a] = b
        apply(study, updates)
        return 0
    key, acc = sys.argv[1], sys.argv[2]
    if key.upper() == "STUDY":
        apply(acc, {})
    else:
        apply(None, {key: acc})
    return 0


def apply(study: str | None, sample_updates: dict[str, str]) -> None:
    rows = []
    with PATH.open() as fh:
        r = csv.DictReader(fh, delimiter="\t")
        fields = r.fieldnames
        for row in r:
            rows.append(row)
    for row in rows:
        if study and (row["library_name"] == "STUDY" or True):
            if study:
                row["study_accession"] = study
        lib = row["library_name"]
        alias = row.get("sample_alias", "")
        if lib in sample_updates:
            row["sample_accession"] = sample_updates[lib]
        if alias in sample_updates:
            row["sample_accession"] = sample_updates[alias]
    with PATH.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=fields, delimiter="\t")
        w.writeheader()
        w.writerows(rows)
    print(f"Updated {PATH}")
    print("Re-run: python3 generate_ena_pack.py")


if __name__ == "__main__":
    sys.exit(main())
