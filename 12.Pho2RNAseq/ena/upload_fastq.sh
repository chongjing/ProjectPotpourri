#!/usr/bin/env bash
# Upload FASTQ files to the Webin FTP dropbox. Does not submit metadata.
# Credentials MUST come from the environment — do not put the password in this file.
#
#   export ENA_WEBIN_USER=Webin-69760
#   export ENA_WEBIN_PASSWORD='...'
#   ./upload_fastq.sh
#
# Files land at: webin-cli/reads/<LIBRARY_NAME>/<FASTQ>
set -euo pipefail

: "${ENA_WEBIN_USER:?set ENA_WEBIN_USER}"
: "${ENA_WEBIN_PASSWORD:?set ENA_WEBIN_PASSWORD}"

HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
HOST="ftp://webin2.ebi.ac.uk"
LOG="$HERE/upload_fastq.log"

if [[ ! -f metadata.tsv ]]; then
  echo "metadata.tsv missing; run generate_ena_pack.py first" >&2
  exit 1
fi

echo "=== upload start $(date -Iseconds) ===" | tee -a "$LOG"

python3 - <<'PY' | while IFS=$'\t' read -r lib fq; do
  remote="${HOST}/webin-cli/reads/${lib}/${fq}"
  echo "=== ${lib}: ${fq} ===" | tee -a "$LOG"
  curl -T "${HERE}/${fq}" \
    --retry 5 --retry-delay 30 \
    --connect-timeout 60 --max-time 3600 \
    --ftp-create-dirs \
    --user "${ENA_WEBIN_USER}:${ENA_WEBIN_PASSWORD}" \
    "$remote" \
    && echo "DONE ${fq}" | tee -a "$LOG" \
    || { echo "FAILED ${fq}" | tee -a "$LOG"; exit 1; }
done
import csv
with open("metadata.tsv") as fh:
    for row in csv.DictReader(fh, delimiter="\t"):
        print(f"{row['library_name']}\t{row['fq1']}")
        print(f"{row['library_name']}\t{row['fq2']}")
PY

echo "=== upload finished $(date -Iseconds) ===" | tee -a "$LOG"
