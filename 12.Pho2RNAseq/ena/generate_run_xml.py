#!/usr/bin/env python3
"""Write experiment/run/submission XML for each library (PAIRED layout).

Used after FASTQs are on webin2.ebi.ac.uk at:
  webin-cli/reads/<LIBRARY_NAME>/<FASTQ>
"""
from __future__ import annotations

import csv
import html
from pathlib import Path

ROOT = Path(__file__).resolve().parent
OUT = ROOT / "xml" / "runs"


def esc(s: str) -> str:
    return html.escape(s, quote=True)


def write_one(row: dict) -> None:
    lib = row["library_name"]
    alias = f"webin-reads-{lib}"
    d = OUT / lib
    d.mkdir(parents=True, exist_ok=True)
    study = row["study_accession"]
    sample = row["sample_accession"]
    title = row["sample_title"]
    remote1 = f"webin-cli/reads/{lib}/{row['fq1']}"
    remote2 = f"webin-cli/reads/{lib}/{row['fq2']}"

    (d / "experiment.xml").write_text(
        f"""<?xml version="1.0" encoding="UTF-8"?>
<EXPERIMENT_SET>
  <EXPERIMENT alias="{esc(alias)}">
    <TITLE>Raw reads: {esc(lib)} — {esc(title)}</TITLE>
    <STUDY_REF accession="{esc(study)}"/>
    <DESIGN>
      <DESIGN_DESCRIPTION>150 bp paired-end mRNA-seq, Illumina NovaSeq 6000, Novogene Cambridge.</DESIGN_DESCRIPTION>
      <SAMPLE_DESCRIPTOR accession="{esc(sample)}"/>
      <LIBRARY_DESCRIPTOR>
        <LIBRARY_NAME>{esc(lib)}</LIBRARY_NAME>
        <LIBRARY_STRATEGY>RNA-Seq</LIBRARY_STRATEGY>
        <LIBRARY_SOURCE>TRANSCRIPTOMIC</LIBRARY_SOURCE>
        <LIBRARY_SELECTION>Oligo-dT</LIBRARY_SELECTION>
        <LIBRARY_LAYOUT>
          <PAIRED/>
        </LIBRARY_LAYOUT>
      </LIBRARY_DESCRIPTOR>
    </DESIGN>
    <PLATFORM>
      <ILLUMINA>
        <INSTRUMENT_MODEL>Illumina NovaSeq 6000</INSTRUMENT_MODEL>
      </ILLUMINA>
    </PLATFORM>
  </EXPERIMENT>
</EXPERIMENT_SET>
"""
    )
    (d / "run.xml").write_text(
        f"""<?xml version="1.0" encoding="UTF-8"?>
<RUN_SET>
  <RUN alias="{esc(alias)}">
    <TITLE>Raw reads: {esc(lib)}</TITLE>
    <EXPERIMENT_REF refname="{esc(alias)}"/>
    <DATA_BLOCK>
      <FILES>
        <FILE filename="{esc(remote1)}" filetype="fastq" checksum_method="MD5" checksum="{row['fq1_md5']}"/>
        <FILE filename="{esc(remote2)}" filetype="fastq" checksum_method="MD5" checksum="{row['fq2_md5']}"/>
      </FILES>
    </DATA_BLOCK>
  </RUN>
</RUN_SET>
"""
    )
    (d / "submission.xml").write_text(
        f"""<?xml version="1.0" encoding="UTF-8"?>
<SUBMISSION alias="submit-reads-{esc(lib)}">
  <ACTIONS>
    <ACTION>
      <ADD/>
    </ACTION>
  </ACTIONS>
</SUBMISSION>
"""
    )


def main() -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    n = 0
    with (ROOT / "metadata.tsv").open() as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if not row["study_accession"] or not row["sample_accession"]:
                raise SystemExit(f"missing accession for {row['library_name']}")
            if "TO_COMPUTE" in row["fq1_md5"]:
                raise SystemExit(f"missing MD5 for {row['library_name']}")
            write_one(row)
            n += 1
    print(f"Wrote {n} experiment/run/submission triples under {OUT}")


if __name__ == "__main__":
    main()
