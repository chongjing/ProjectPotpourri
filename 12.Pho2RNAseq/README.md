# Rice AM × PHO2 root RNA-seq (Saskia)

Rice (*Oryza sativa* ssp. Japonica cv. Nipponbare) root mRNA-seq: **wild type vs `pho2`**, **low vs high phosphate**, **mock vs *Rhizophagus irregularis***, four biological replicates (32 libraries).

Paper title: *PHOSPHATE OVERACCUMULATOR 2 (PHO2) is a negative regulator of arbuscular mycorrhizal symbiosis* (`docs/Project.txt`).

This GitHub folder is the **curated analysis + ENA submission pack**. Raw FASTQ and BAMs stay on CSD3 (`/rds/user/cx264/hpc-work/project/2.Jeongmin/01.Saskia`) and at ENA.

---

## Experimental design

| Libraries | Genotype | Pi | Inoculation |
|-----------|----------|----|-------------|
| r2310_1–4 | WT | 25 µM (LP) | mock |
| r2310_5–8 | WT | 25 µM (LP) | *R. irregularis* (600 spores) |
| r2310_9–12 | pho2 | 25 µM (LP) | mock |
| r2310_13–16 | pho2 | 25 µM (LP) | *R. irregularis* |
| r2310_17–20 | WT | 250 µM (HP) | mock |
| r2310_21–24 | WT | 250 µM (HP) | *R. irregularis* |
| r2310_25–28 | pho2 | 250 µM (HP) | mock |
| r2310_29–32 | pho2 | 250 µM (HP) | *R. irregularis* |

Growth: autoclaved silica sand, CONVIRON EVO (12 h light, 28/20 °C, 65% RH), half-strength Hoagland. Tissue: **root**. Collection: **2023-10**, United Kingdom. Sequenced 150 bp PE mRNA (oligo-dT) on **Illumina NovaSeq 6000** (Novogene, Cambridge). `pho2` allele not distinguished (`pho2-1` NF2586 vs `pho2-2` NE9017).

Concatenated two-flowcell libraries (one ENA run each): **r2310_1, 8, 19, 23**.

---

## ENA accessions (PRIVATE until 2027-09-05)

| Object | Accession |
|--------|-----------|
| BioProject | [PRJEB125245](https://www.ebi.ac.uk/ena/browser/view/PRJEB125245) |
| Study | ERP204980 |
| Samples | ERS31241328–ERS31241359 |
| Experiments | ERX17188148–ERX17188179 |
| Runs | ERR17796168–ERR17796199 |

Full table: `ena/accessions_registered.tsv`. Paper text: `ena/data_availability_statement.md`. How the submission was done: `ena/PIPELINE.md`.

---

## Files

| Path | What |
|------|------|
| `docs/Project.txt` | Title and abstract |
| `docs/MaterialMethod.txt` | Methods (growth, RNA-seq, DE) |
| `docs/11.script.txt` | HPC lab notebook (historical; includes later Arabidopsis notes) |
| `ena/` | Webin XML, manifests, upload/submit scripts, checksums |
| `de/FDR_0.01/` | edgeR QL F-tests: **FDR &lt; 0.05 and \|logFC\| &gt; 1** (analysis of record) |
| `de/FDR_0.05/` | Same 12 contrasts called on unadjusted *P* &lt; 0.05 (exploratory lists only) |
| `de/markers/1.gene.list` | AM / Pi marker genes (PHO2, PT11, AM1/3/14, IPS1, SQD2.1, PAP10) |

HPC mapping used Novoalign + HTSeq (`--type transcript`). The paper methods use Salmon + DESeq2. Raw reads are the same.

---

## Differential expression (12 pairwise tests)

edgeR `glmQLFit` / `glmQLFTest`, design `~Time + group` (Time = paired experimental batch 1–4).

| Folder | Contrast |
|--------|----------|
| 01 | WT LP myc vs mock |
| 02 | pho2 LP myc vs mock |
| 03 | WT HP myc vs mock |
| 04 | pho2 HP myc vs mock |
| 05 | WT mock HP vs LP |
| 06 | pho2 mock HP vs LP |
| 07 | WT myc HP vs LP |
| 08 | pho2 myc HP vs LP |
| 09 | LP mock WT vs pho2 |
| 10 | LP myc WT vs pho2 |
| 11 | HP mock WT vs pho2 |
| 12 | HP myc WT vs pho2 |

---

## Reproduce ENA submit (credentials not in this repo)

```bash
export ENA_WEBIN_USER=Webin-69760
export ENA_WEBIN_PASSWORD='...'   # never commit
cd ena
python3 generate_ena_pack.py
# FASTQs are not in git; point metadata paths at HPC or re-download from ENA
```
