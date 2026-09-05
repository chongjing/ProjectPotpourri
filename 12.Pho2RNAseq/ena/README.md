# ENA submission pack — Saskia rice AM × PHO2 RNA-seq

Local files only. Nothing has been uploaded or registered at ENA yet.

| | |
|---|---|
| Study alias | `Saskia_PHO2_AM_RNAseq` |
| BioProject | **PRJEB125245** (study ERP204980) |
| Samples | **ERS31241328–ERS31241359** (BioSamples SAMEA123420480–SAMEA123420511) |
| Experiments / runs | **ERX17188148–ERX17188179** / **ERR17796168–ERR17796199** (32 paired FASTQ runs, PRIVATE) |
| Hold date | **2027-09-05** (PRIVATE) |
| Webin account | Webin-69760 |
| Samples | 32 (one BioSample + one paired FASTQ run each) |
| Data | 64 gzip FASTQ, ~146 GB, **symlinked** (not copied) |
| Instrument | Illumina NovaSeq 6000 |
| Library | mRNA oligo-dT, 150 bp paired-end, Novogene Cambridge |

Do **not** reuse `PRJEB114114` (that is a different organism/study).

---

## What is already done here

1. `generate_ena_pack.py` writes metadata, plant-sample TSV, webin-cli manifests, checksums, and FASTQ symlinks.
2. Four libraries (`r2310_1`, `r2310_8`, `r2310_19`, `r2310_23`) were concatenated from two flow cells. Their MD5s are **recomputed** from the files on disk. The other 28 use Novogene `MD5.txt` where the filename still matches.
3. Manifests currently use **aliases** (`Saskia_PHO2_AM_RNAseq`, `Saskia_r2310_N`). After you register study/samples in Webin with those exact aliases, webin-cli can resolve them. You can later swap in `PRJEB`/`ERS` via `accessions.tsv`.

---

## Registered (2026-09-05)

Study and 32 plant samples were registered on the **production** Webin drop-box (checklist ERC000037). Status PRIVATE until 2027-09-05. Full table: `accessions_registered.tsv`. Receipts: `xml/receipt_prod_study.xml`, `xml/receipt_prod_sample.xml`.

The public browser page will be empty until release:
https://www.ebi.ac.uk/ena/browser/view/PRJEB125245

Read manifests now contain `STUDY PRJEB125245` and the matching `SAMPLE ERS…`.

### 3. FASTQ upload and ERR/ERX

Full procedure: **`PIPELINE.md`**.

Slurm upload job: `upload_fastq.srun` (job id in `upload_job.id`). After 64 files are on webin2, `submit_runs_rest.sh` registers PAIRED experiments/runs.

If webin-cli FTP times out (it did on the previous 8–9 GB Sebastian job), `upload_fastq.sh` is the same curl workaround. Files here are 3–6 GB so native webin-cli may work; still run curl first if in doubt.

If webin-cli writes `<SINGLE />` for paired data, stop and generate experiment/run XML as in the Sebastian pack (`ena_xmls/`). Do not accept SINGLE layout.

---

## Sample map

| Libraries | Genotype | Pi | Inoculation |
|---|---|---|---|
| r2310_1–4 | WT | 25 uM LP | mock |
| r2310_5–8 | WT | 25 uM LP | *R. irregularis* |
| r2310_9–12 | pho2 | 25 uM LP | mock |
| r2310_13–16 | pho2 | 25 uM LP | *R. irregularis* |
| r2310_17–20 | WT | 250 uM HP | mock |
| r2310_21–24 | WT | 250 uM HP | *R. irregularis* |
| r2310_25–28 | pho2 | 250 uM HP | mock |
| r2310_29–32 | pho2 | 250 uM HP | *R. irregularis* |

Concatenated two-flowcell libraries (one run each): **r2310_1, 8, 19, 23**.

The original unmerged lane FASTQs still exist as **hidden** files in the sample folders (names starting with `.`). They are not submitted. If you later want two runs per library, those files can be un-hidden and registered as extra runs on the same BioSample.

---

## Files

| File | Role |
|---|---|
| `generate_ena_pack.py` | Rebuild metadata, manifests, symlinks |
| `metadata.tsv` | Master table |
| `sample_checklist_ERC000037.tsv` | Webin sample spreadsheet source |
| `study_text.md` | Title/abstract/hold date for the portal |
| `accessions.tsv` | Fill PRJEB/ERS after registration |
| `manifests/` | webin-cli `-context reads` |
| `checksums.md5` | MD5 of the 64 FASTQ files |
| `upload_fastq.sh` | curl FTP (env credentials) |
| `submit_reads.sh` | webin-cli validate/submit (env credentials) |
| `fill_accessions.py` | Write PRJEB/ERS into accessions.tsv |
| `*.fq.gz` | Symlinks to `../X204SC…/01.RawData/` |

No Webin password is stored in this directory.
