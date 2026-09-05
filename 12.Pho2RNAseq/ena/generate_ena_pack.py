#!/usr/bin/env python3
"""Build the local ENA submission pack for Saskia rice AM x PHO2 RNA-seq.

Does not upload or talk to ENA. STUDY/SAMPLE in manifests are Webin *aliases*
that must be used unchanged when registering the study and samples in the
Webin portal.
"""
from __future__ import annotations

import csv
import hashlib
import os
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parent
DATA = ROOT.parent
PACK01 = DATA / "X204SC23091440-Z01-F002_01" / "01.RawData"
PACK02 = DATA / "X204SC23091440-Z01-F002_02" / "01.RawData"

STUDY_ALIAS = "Saskia_PHO2_AM_RNAseq"
HOLD_DATE = "2027-09-05"
INSTRUMENT = "Illumina NovaSeq 6000"
LIBRARY_SOURCE = "TRANSCRIPTOMIC"
LIBRARY_SELECTION = "Oligo-dT"
LIBRARY_STRATEGY = "RNA-Seq"

# Concatenated on 2024-01-23 from two NovaSeq flow cells. Novogene MD5.txt
# still lists the original four lane files, which are no longer on disk.
CONCATENATED = {1, 8, 19, 23}

# Novogene delivery split.
IN_PACK01 = {4, 5, 6, 7, 8, 9, 32}

CONDITIONS = {
    1: ("WT", "LP", "mock", 1),
    2: ("WT", "LP", "mock", 2),
    3: ("WT", "LP", "mock", 3),
    4: ("WT", "LP", "mock", 4),
    5: ("WT", "LP", "myc", 1),
    6: ("WT", "LP", "myc", 2),
    7: ("WT", "LP", "myc", 3),
    8: ("WT", "LP", "myc", 4),
    9: ("pho2", "LP", "mock", 1),
    10: ("pho2", "LP", "mock", 2),
    11: ("pho2", "LP", "mock", 3),
    12: ("pho2", "LP", "mock", 4),
    13: ("pho2", "LP", "myc", 1),
    14: ("pho2", "LP", "myc", 2),
    15: ("pho2", "LP", "myc", 3),
    16: ("pho2", "LP", "myc", 4),
    17: ("WT", "HP", "mock", 1),
    18: ("WT", "HP", "mock", 2),
    19: ("WT", "HP", "mock", 3),
    20: ("WT", "HP", "mock", 4),
    21: ("WT", "HP", "myc", 1),
    22: ("WT", "HP", "myc", 2),
    23: ("WT", "HP", "myc", 3),
    24: ("WT", "HP", "myc", 4),
    25: ("pho2", "HP", "mock", 1),
    26: ("pho2", "HP", "mock", 2),
    27: ("pho2", "HP", "mock", 3),
    28: ("pho2", "HP", "mock", 4),
    29: ("pho2", "HP", "myc", 1),
    30: ("pho2", "HP", "myc", 2),
    31: ("pho2", "HP", "myc", 3),
    32: ("pho2", "HP", "myc", 4),
}

PI_UM = {"LP": "25", "HP": "250"}
PI_LABEL = {"LP": "low phosphate (25 uM)", "HP": "high phosphate (250 uM)"}
INOC_LABEL = {
    "mock": "mock inoculation (water)",
    "myc": "Rhizophagus irregularis (600 spores)",
}
GENO_LABEL = {
    "WT": "wild type",
    "pho2": "pho2 Tos17 insertion mutant of OsPHO2 (LOC_Os05g48390); allele not distinguished between pho2-1 (NF2586) and pho2-2 (NE9017)",
}


def sample_dir(n: int) -> Path:
    pack = PACK01 if n in IN_PACK01 else PACK02
    d = pack / f"r2310_{n}"
    if not d.is_dir():
        raise FileNotFoundError(d)
    return d


def fastq_pair(n: int) -> tuple[Path, Path]:
    d = sample_dir(n)
    # Ignore leftover hidden lane files (e.g. .r2310_1_...HWT7HDSX7_L2_1.fq.gz).
    # The submitted objects are the visible concatenated/delivered FASTQs.
    fqs = sorted(p for p in d.glob("*.fq.gz") if not p.name.startswith("."))
    r1 = [p for p in fqs if p.name.endswith("_1.fq.gz")]
    r2 = [p for p in fqs if p.name.endswith("_2.fq.gz")]
    if len(r1) != 1 or len(r2) != 1:
        raise RuntimeError(f"r2310_{n}: expected one R1 and one R2, found {fqs}")
    return r1[0], r2[0]


def novogene_md5(n: int) -> dict[str, str]:
    p = sample_dir(n) / "MD5.txt"
    out = {}
    for line in p.read_text().splitlines():
        line = line.strip()
        if not line:
            continue
        md5, name = line.split()
        out[name] = md5
    return out


def md5_file(path: Path, chunk: int = 1024 * 1024) -> str:
    h = hashlib.md5()
    with path.open("rb") as fh:
        while True:
            b = fh.read(chunk)
            if not b:
                break
            h.update(b)
    return h.hexdigest()


def title(n: int) -> str:
    geno, pi, inoc, rep = CONDITIONS[n]
    g = "WT" if geno == "WT" else "pho2"
    return (
        f"Oryza sativa Nipponbare {g} root, {PI_LABEL[pi]}, "
        f"{INOC_LABEL[inoc]}, replicate {rep}"
    )


def description(n: int) -> str:
    geno, pi, inoc, rep = CONDITIONS[n]
    bits = [
        "Rice (Oryza sativa L. ssp. Japonica cv. Nipponbare) root RNA-seq.",
        f"Genotype: {GENO_LABEL[geno]}.",
        "Grown in autoclaved silica 14/25 sand in a CONVIRON EVO chamber "
        "(12 h light/12 h dark, 350 uE m-2 s-1, 28/20 C, 65% RH).",
        "Seedlings watered with RO water for one week, then fertilised twice "
        f"weekly with half-strength Hoagland solution at {PI_UM[pi]} uM phosphate "
        f"({pi}) plus 0.01% (w/v) Sequestrene Rapid iron supplement.",
        (
            "Mock inoculated with RO water."
            if inoc == "mock"
            else "Inoculated with 600 freshly extracted Rhizophagus irregularis spores from hairy carrot root organ cultures."
        ),
        f"Biological replicate {rep} (experimental batch {rep}).",
        "Organism part: root. Collection date: 2023-10. Geographic location: United Kingdom "
        "(Crop Science Centre, University of Cambridge).",
        "mRNA library (oligo-dT), 150 bp paired-end, Illumina NovaSeq 6000, Novogene (Cambridge, UK).",
        f"Library name: r2310_{n}.",
    ]
    if n in CONCATENATED:
        bits.append(
            "This library was sequenced on two NovaSeq 6000 flow cells "
            "(HWT7HDSX7 and H2WFLDSXC); the two lanes were concatenated into "
            "a single FASTQ pair before submission."
        )
    return " ".join(bits)


def sample_alias(n: int) -> str:
    return f"Saskia_r2310_{n}"


def write_metadata(rows: list[dict]) -> None:
    cols = [
        "library_name",
        "sample_alias",
        "study_alias",
        "genotype",
        "phosphate",
        "phosphate_uM",
        "inoculation",
        "replicate",
        "concatenated_lanes",
        "instrument",
        "library_source",
        "library_selection",
        "library_strategy",
        "layout",
        "fq1",
        "fq1_md5",
        "fq1_bytes",
        "fq2",
        "fq2_md5",
        "fq2_bytes",
        "fq1_src",
        "fq2_src",
        "sample_title",
        "sample_description",
        "md5_source",
        "study_accession",
        "sample_accession",
    ]
    path = ROOT / "metadata.tsv"
    with path.open("w", newline="") as fh:
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t", extrasaction="ignore")
        w.writeheader()
        for row in rows:
            w.writerow(row)
    print(f"Wrote {path}")


def write_sample_checklist(rows: list[dict]) -> None:
    """TSV to paste into the Webin ERC000037 plant sample spreadsheet."""
    path = ROOT / "sample_checklist_ERC000037.tsv"
    # Columns that exist on ERC000037 plus the three Webin identity fields.
    # Phosphate / inoculation / replicate are in sample_title and sample_description
    # so they are not extra spreadsheet columns (Webin rejects unknown headers).
    cols = [
        "sample_alias",
        "tax_id",
        "scientific_name",
        "common_name",
        "sample_title",
        "sample_description",
        "collection date",
        "geographic location (country and/or sea)",
        "geographic location (region and locality)",
        "cultivar",
        "genotype",
        "organism part",
        "isolation_source",
        "infect",
        "collected_by",
    ]
    with path.open("w", newline="") as fh:
        fh.write("#checklist_accession\tERC000037\n")
        fh.write(
            "# Paste these rows into the Webin Portal 'Submit Samples' spreadsheet "
            "after selecting ENA Plant Sample Checklist (ERC000037). "
            "Keep sample_alias exactly as written. Hold date is set on the STUDY, not here.\n"
        )
        w = csv.DictWriter(fh, fieldnames=cols, delimiter="\t")
        w.writeheader()
        for row in rows:
            n = int(row["library_name"].split("_")[1])
            geno, pi, inoc, rep = CONDITIONS[n]
            w.writerow(
                {
                    "sample_alias": row["sample_alias"],
                    "tax_id": "39947",
                    "scientific_name": "Oryza sativa Japonica Group",
                    "common_name": "rice",
                    "sample_title": row["sample_title"],
                    "sample_description": row["sample_description"],
                    "collection date": "2023-10",
                    "geographic location (country and/or sea)": "United Kingdom",
                    "geographic location (region and locality)": "Cambridge",
                    "cultivar": "Nipponbare",
                    "genotype": "wild type" if geno == "WT" else "pho2",
                    "organism part": "root",
                    "isolation_source": "root",
                    "infect": "none" if inoc == "mock" else "Rhizophagus irregularis",
                    "collected_by": "Jeongmin Choi",
                }
            )
    print(f"Wrote {path}")


def write_manifest(row: dict) -> None:
    man_dir = ROOT / "manifests"
    man_dir.mkdir(exist_ok=True)
    name = row["library_name"]
    study = row["study_accession"] or STUDY_ALIAS
    sample = row["sample_accession"] or row["sample_alias"]
    text = (
        f"STUDY\t{study}\n"
        f"SAMPLE\t{sample}\n"
        f"NAME\t{name}\n"
        f"INSTRUMENT\t{INSTRUMENT}\n"
        f"LIBRARY_NAME\t{name}\n"
        f"LIBRARY_SOURCE\t{LIBRARY_SOURCE}\n"
        f"LIBRARY_SELECTION\t{LIBRARY_SELECTION}\n"
        f"LIBRARY_STRATEGY\t{LIBRARY_STRATEGY}\n"
        f"DESCRIPTION\t{row['sample_title']}; 150 bp paired-end mRNA, Novogene Cambridge.\n"
        f"FASTQ\t{row['fq1']}\n"
        f"FASTQ\t{row['fq2']}\n"
    )
    path = man_dir / f"{name}.manifest.txt"
    path.write_text(text)
    print(f"Wrote {path}")


def write_checksums(rows: list[dict]) -> None:
    lines = []
    for row in rows:
        lines.append(f"{row['fq1_md5']}  {row['fq1']}")
        lines.append(f"{row['fq2_md5']}  {row['fq2']}")
    path = ROOT / "checksums.md5"
    path.write_text("\n".join(lines) + "\n")
    print(f"Wrote {path}")


def write_accessions_template() -> None:
    path = ROOT / "accessions.tsv"
    if path.exists():
        print(f"Keeping existing {path}")
        return
    with path.open("w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["library_name", "sample_alias", "sample_accession", "study_accession", "run_accession", "experiment_accession"])
        w.writerow(["STUDY", STUDY_ALIAS, "", "PRJEB_TO_FILL", "", ""])
        for n in range(1, 33):
            w.writerow([f"r2310_{n}", sample_alias(n), "ERS_TO_FILL", "PRJEB_TO_FILL", "", ""])
    print(f"Wrote {path}")


def load_accessions() -> tuple[str, dict[str, str]]:
    path = ROOT / "accessions.tsv"
    study = ""
    sample_acc = {}
    if not path.exists():
        return study, sample_acc
    with path.open() as fh:
        r = csv.DictReader(fh, delimiter="\t")
        for row in r:
            acc = (row.get("sample_accession") or "").strip()
            study_acc = (row.get("study_accession") or "").strip()
            if row["library_name"] == "STUDY":
                if study_acc and not study_acc.startswith("PRJEB_TO_FILL"):
                    study = study_acc
                continue
            if acc and not acc.startswith("ERS_TO_FILL"):
                sample_acc[row["library_name"]] = acc
            if study_acc and not study_acc.startswith("PRJEB_TO_FILL"):
                study = study_acc
    return study, sample_acc


def ensure_symlink(src: Path, dst: Path) -> None:
    if dst.is_symlink() or dst.exists():
        if dst.is_symlink() and dst.resolve() == src.resolve():
            return
        dst.unlink()
    os.symlink(src, dst)


def main() -> int:
    compute_missing = "--compute-md5" in sys.argv
    study_acc, sample_acc = load_accessions()
    rows = []
    missing_md5 = []

    for n in range(1, 33):
        geno, pi, inoc, rep = CONDITIONS[n]
        r1, r2 = fastq_pair(n)
        known = novogene_md5(n)
        md5_src = "novogene"
        m1 = known.get(r1.name)
        m2 = known.get(r2.name)
        if not m1 or not m2:
            # Reuse previously computed MD5s from metadata.tsv if present.
            old = {}
            oldp = ROOT / "metadata.tsv"
            if oldp.exists():
                with oldp.open() as fh:
                    for prev in csv.DictReader(fh, delimiter="\t"):
                        if prev["library_name"] == f"r2310_{n}":
                            old = prev
                            break
            if (
                old.get("fq1") == r1.name
                and old.get("fq2") == r2.name
                and old.get("fq1_md5")
                and old.get("fq2_md5")
                and "TO_COMPUTE" not in old["fq1_md5"]
            ):
                m1, m2 = old["fq1_md5"], old["fq2_md5"]
                md5_src = old.get("md5_source") or "computed"
            elif compute_missing:
                md5_src = "computed"
                print(f"Computing MD5 r2310_{n} (concatenated or renamed) ...", flush=True)
                m1 = md5_file(r1)
                m2 = md5_file(r2)
            else:
                md5_src = "MISSING"
                m1 = m1 or "MD5_TO_COMPUTE"
                m2 = m2 or "MD5_TO_COMPUTE"
                missing_md5.append(n)
        dst1 = ROOT / r1.name
        dst2 = ROOT / r2.name
        ensure_symlink(r1, dst1)
        ensure_symlink(r2, dst2)
        lib = f"r2310_{n}"
        rows.append(
            {
                "library_name": lib,
                "sample_alias": sample_alias(n),
                "study_alias": STUDY_ALIAS,
                "genotype": geno,
                "phosphate": pi,
                "phosphate_uM": PI_UM[pi],
                "inoculation": inoc,
                "replicate": str(rep),
                "concatenated_lanes": "yes" if n in CONCATENATED else "no",
                "instrument": INSTRUMENT,
                "library_source": LIBRARY_SOURCE,
                "library_selection": LIBRARY_SELECTION,
                "library_strategy": LIBRARY_STRATEGY,
                "layout": "PAIRED",
                "fq1": r1.name,
                "fq1_md5": m1,
                "fq1_bytes": str(r1.stat().st_size),
                "fq2": r2.name,
                "fq2_md5": m2,
                "fq2_bytes": str(r2.stat().st_size),
                "fq1_src": str(r1),
                "fq2_src": str(r2),
                "sample_title": title(n),
                "sample_description": description(n),
                "md5_source": md5_src,
                "study_accession": study_acc,
                "sample_accession": sample_acc.get(lib, ""),
            }
        )

    write_metadata(rows)
    write_sample_checklist(rows)
    write_checksums(rows)
    write_accessions_template()
    for row in rows:
        write_manifest(row)

    n_link = len(list(ROOT.glob("*.fq.gz")))
    print(f"Symlinks in {ROOT}: {n_link}")
    if missing_md5:
        print(
            "MD5 missing for concatenated/renamed libraries: "
            + ", ".join(f"r2310_{n}" for n in missing_md5)
        )
        print("Re-run: python3 generate_ena_pack.py --compute-md5")
        return 2
    return 0


if __name__ == "__main__":
    sys.exit(main())
