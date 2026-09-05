#!/usr/bin/env python3
"""Write project.xml, sample.xml, submission.xml for Webin drop-box registration.

Study TITLE/DESCRIPTION are capped at 249 characters (ENA INSDC limit).
Sample checklist: ERC000037 (plant). Extra attributes are allowed.
"""
from __future__ import annotations

import csv
import html
from pathlib import Path

ROOT = Path(__file__).resolve().parent
XML = ROOT / "xml"
STUDY_ALIAS = "Saskia_PHO2_AM_RNAseq"
HOLD_DATE = "2027-09-05"

TITLE = (
    "PHOSPHATE OVERACCUMULATOR 2 (PHO2) is a negative regulator "
    "of arbuscular mycorrhizal symbiosis"
)
# ENA requires TITLE and DESCRIPTION length 21-249.
DESCRIPTION = (
    "Root mRNA-seq of Oryza sativa Nipponbare WT and pho2 at 25 or 250 uM "
    "phosphate, mock or Rhizophagus irregularis; 32 libraries, 4 replicates. "
    "Crop Science Centre, University of Cambridge."
)

GROWTH = (
    "Seeds sterilised, germinated 3-6 d at 30 C, grown in autoclaved silica "
    "14/25 sand in a CONVIRON EVO chamber (12 h light, 28/20 C, 65% RH) with "
    "half-strength Hoagland solution at 25 or 250 uM phosphate."
)


def esc(s: str) -> str:
    return html.escape(s, quote=True)


def attr(tag: str, value: str, units: str | None = None) -> str:
    u = f"\n        <UNITS>{esc(units)}</UNITS>" if units else ""
    return (
        f"      <SAMPLE_ATTRIBUTE>\n"
        f"        <TAG>{esc(tag)}</TAG>\n"
        f"        <VALUE>{esc(value)}</VALUE>{u}\n"
        f"      </SAMPLE_ATTRIBUTE>"
    )


def write_project() -> None:
    assert 20 < len(TITLE) < 250, len(TITLE)
    assert 20 < len(DESCRIPTION) < 250, len(DESCRIPTION)
    xml = f"""<?xml version="1.0" encoding="UTF-8"?>
<PROJECT_SET>
  <PROJECT alias="{esc(STUDY_ALIAS)}">
    <NAME>{esc(STUDY_ALIAS)}</NAME>
    <TITLE>{esc(TITLE)}</TITLE>
    <DESCRIPTION>{esc(DESCRIPTION)}</DESCRIPTION>
    <SUBMISSION_PROJECT>
      <SEQUENCING_PROJECT/>
    </SUBMISSION_PROJECT>
  </PROJECT>
</PROJECT_SET>
"""
    (XML / "project.xml").write_text(xml)
    print(f"Wrote {XML / 'project.xml'} title={len(TITLE)} desc={len(DESCRIPTION)}")


def write_submission(name: str, hold: bool) -> None:
    hold_xml = ""
    if hold:
        hold_xml = f"""
    <ACTION>
      <HOLD HoldUntilDate="{HOLD_DATE}"/>
    </ACTION>"""
    xml = f"""<?xml version="1.0" encoding="UTF-8"?>
<SUBMISSION alias="{esc(name)}">
  <CONTACTS>
    <CONTACT name="Jeongmin Choi" inform_on_error="jc913@cam.ac.uk" inform_on_status="jc913@cam.ac.uk"/>
    <CONTACT name="Chongjing Xia" inform_on_error="xiachongjing@gmail.com" inform_on_status="xiachongjing@gmail.com"/>
  </CONTACTS>
  <ACTIONS>
    <ACTION>
      <ADD/>
    </ACTION>{hold_xml}
  </ACTIONS>
</SUBMISSION>
"""
    (XML / f"{name}.xml").write_text(xml)
    print(f"Wrote {XML / f'{name}.xml'}")


def write_samples() -> None:
    rows = list(csv.DictReader((ROOT / "metadata.tsv").open(), delimiter="\t"))
    blocks = []
    for row in rows:
        geno = "wild type" if row["genotype"] == "WT" else "pho2"
        infect = (
            "none"
            if row["inoculation"] == "mock"
            else "Rhizophagus irregularis"
        )
        biotic = "free living" if row["inoculation"] == "mock" else "symbiont"
        attrs = [
            attr("ENA-CHECKLIST", "ERC000037"),
            attr("collection date", "2023-10"),
            attr("geographic location (country and/or sea)", "United Kingdom"),
            attr("geographic location (region and locality)", "Cambridge"),
            attr("geographic location (latitude)", "not collected", "DD"),
            attr("geographic location (longitude)", "not collected", "DD"),
            attr("plant structure", "root"),
            attr("plant developmental stage", "vegetative"),
            attr("plant growth medium", "silica sand"),
            attr("isolation and growth condition", GROWTH),
            attr("growth facility", "growth chamber"),
            attr("organism common name", "rice"),
            attr("genotype", geno),
            attr("cultivar", "Nipponbare"),
            attr("collected_by", "Jeongmin Choi"),
            attr("sample health state", "healthy"),
            attr("observed biotic relationship", biotic),
            attr("infect", infect),
            attr("phosphate treatment", f"{row['phosphate_uM']} uM {row['phosphate']}"),
            attr("inoculation", row["inoculation"]),
            attr("replicate", row["replicate"]),
        ]
        desc = row["sample_description"]
        blocks.append(
            f"""  <SAMPLE alias="{esc(row['sample_alias'])}">
    <TITLE>{esc(row['sample_title'])}</TITLE>
    <SAMPLE_NAME>
      <TAXON_ID>39947</TAXON_ID>
      <SCIENTIFIC_NAME>Oryza sativa Japonica Group</SCIENTIFIC_NAME>
      <COMMON_NAME>rice</COMMON_NAME>
    </SAMPLE_NAME>
    <DESCRIPTION>{esc(desc)}</DESCRIPTION>
    <SAMPLE_ATTRIBUTES>
{chr(10).join(attrs)}
    </SAMPLE_ATTRIBUTES>
  </SAMPLE>"""
        )
    xml = (
        '<?xml version="1.0" encoding="UTF-8"?>\n'
        "<SAMPLE_SET>\n"
        + "\n".join(blocks)
        + "\n</SAMPLE_SET>\n"
    )
    (XML / "sample.xml").write_text(xml)
    print(f"Wrote {XML / 'sample.xml'} n={len(blocks)}")
    # one-sample subset for the test server
    first = blocks[0]
    (XML / "sample_test.xml").write_text(
        '<?xml version="1.0" encoding="UTF-8"?>\n<SAMPLE_SET>\n'
        + first
        + "\n</SAMPLE_SET>\n"
    )


def main() -> None:
    XML.mkdir(exist_ok=True)
    write_project()
    write_samples()
    write_submission("submission_study", hold=True)
    write_submission("submission_samples", hold=True)


if __name__ == "__main__":
    main()
