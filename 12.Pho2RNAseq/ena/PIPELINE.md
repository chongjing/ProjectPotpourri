# ENA raw-read pipeline (Saskia rice AM × PHO2)

This is the reproducible path used on CSD3. Do **not** copy FASTQs; they are symlinked.

```
FASTQ (symlinks)
    → curl FTP  →  webin2.ebi.ac.uk/webin-cli/reads/<lib>/<fq>
    → REST XML  →  EXPERIMENT (ERX) + RUN (ERR)
    → study PRJEB125245  (HOLD 2027-09-05)
```

Why not webin-cli upload: on this HPC, webin-cli 9.0.3’s built-in FTP timed out on multi-GB files (Sebastian pack). Files here are 3–6 GB; still use curl. webin-cli is kept for optional `-validate` only.

---

## Status (2026-09-05)

All steps complete. Data are **PRIVATE until 2027-09-05**.

| Object | Accession |
|---|---|
| BioProject | **PRJEB125245** (ERP204980) |
| Samples | ERS31241328–ERS31241359 |
| Experiments | ERX17188148–ERX17188179 |
| Runs | **ERR17796168–ERR17796199** |
| FASTQ upload | Slurm 34923681, 64/64 files, 1 h 15 min |

Full table: `accessions_registered.tsv`. Paper text: `data_availability_statement.md`.

---

## 0. Already done

| Step | Result |
|---|---|
| Study registered | **PRJEB125245** / ERP204980 |
| 32 samples ERC000037 | **ERS31241328–ERS31241359** |
| Hold | 2027-09-05 PRIVATE |
| Manifests | `manifests/r2310_*.manifest.txt` (STUDY/SAMPLE filled) |
| MD5 | Novogene for 28 libraries; recomputed for concatenated r2310_1, 8, 19, 23 |

---

## 1. Upload FASTQs (~146 GB, 64 files)

```bash
cd /rds/user/cx264/hpc-work/project/2.Jeongmin/01.Saskia/1.data/ENA_submit
sbatch upload_fastq.srun
# or, with env credentials:
#   export ENA_WEBIN_USER=Webin-69760
#   export ENA_WEBIN_PASSWORD='...'
#   python3 upload_fastq.py
```

- Resumable: `upload_state.tsv` skips files already `done`.
- 3 attempts × curl `--retry 5`, `--max-time 10800` per file.
- Log: `upload_fastq.log` and Slurm `upload_slurm.o` / `upload_slurm.e`.
- Current Slurm job id: see `upload_job.id`.
- Remote path: `webin-cli/reads/<library_name>/<fastq>` — this is the path in `run.xml`.

Expect ~1.5–6 h at 20–30 MB/s. Re-run the same command to resume.

---

## 2. Register experiments and runs (ERR / ERX)

After `upload_state.tsv` has 64 `done` rows:

```bash
export ENA_WEBIN_USER=Webin-69760
export ENA_WEBIN_PASSWORD='...'
python3 generate_run_xml.py          # xml/runs/r2310_*/{experiment,run,submission}.xml
./submit_runs_rest.sh                # drop-box REST, PAIRED layout
python3 parse_run_receipts.py        # accessions.tsv + data_availability_statement.md
```

Each library is one experiment + one run (two FASTQs). Concatenated libraries are still one run.

Layout is `<PAIRED/>` (webin-cli previously wrote `<SINGLE/>` on this cluster).

Receipts: `xml/runs/r2310_*/receipt.xml` and `xml/run_receipts.tsv`.

---

## 3. Optional: webin-cli validate

Does not replace REST submit.

```bash
export ENA_WEBIN_USER=Webin-69760
export ENA_WEBIN_PASSWORD='...'
./submit_reads.sh validate
```

---

## Scripts

| File | Role |
|---|---|
| `generate_ena_pack.py` | metadata, manifests, symlinks, checksums |
| `generate_registration_xml.py` | study + sample XML (already submitted) |
| `register_study_samples.sh` | drop-box ADD+HOLD for project/samples |
| `upload_fastq.py` | resumable FTP |
| `upload_fastq.srun` | Slurm wrapper (24 h, sapphire) |
| `generate_run_xml.py` | PAIRED experiment/run XML |
| `submit_runs_rest.sh` | REST submit runs |
| `parse_run_receipts.py` | merge ERR/ERX, write DAS |
| `submit_reads.sh` | webin-cli validate/submit (legacy) |

Credentials are never stored in this directory. Export `ENA_WEBIN_USER` and `ENA_WEBIN_PASSWORD` before `sbatch` or REST submit.

---

## Restart / failure

- Upload died: `sbatch upload_fastq.srun` again (skips `done`).
- One REST run failed: inspect `xml/runs/<lib>/receipt.xml`, fix XML, re-run curl for that library only (change alias if ENA already consumed it).
- Wrong layout: do not use webin-cli `-submit`; REST XML is the source of truth.
