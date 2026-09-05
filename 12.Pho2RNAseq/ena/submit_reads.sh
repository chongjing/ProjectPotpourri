#!/usr/bin/env bash
# Validate, then submit, each webin-cli reads manifest.
# STUDY and SAMPLE in manifests must already exist in Webin (alias or accession).
# FASTQ files must already be in this directory (symlinks) and, if using the
# REST/FTP workaround, on webin2.ebi.ac.uk under webin-cli/reads/<NAME>/.
#
#   export ENA_WEBIN_USER=Webin-69760
#   export ENA_WEBIN_PASSWORD='...'
#   ./submit_reads.sh            # validate + submit
#   ./submit_reads.sh validate   # validate only
set -euo pipefail

: "${ENA_WEBIN_USER:?set ENA_WEBIN_USER}"
: "${ENA_WEBIN_PASSWORD:?set ENA_WEBIN_PASSWORD}"

MODE="${1:-submit}"
HERE="$(cd "$(dirname "$0")" && pwd)"
cd "$HERE"
JAVA="/rds/project/rds-FTKWLWDeHys/programs/java/jdk-17.0.10/bin/java"
JAR="$HOME/rds/rds-csc_programmes-FTKWLWDeHys/programs/webin-cli/webin-cli-9.0.3.jar"
OUTDIR="$HERE/webin_out"
mkdir -p "$OUTDIR"

if [[ ! -x "$JAVA" ]]; then
  echo "Java not found: $JAVA" >&2
  exit 1
fi
if [[ ! -f "$JAR" ]]; then
  echo "webin-cli jar not found: $JAR" >&2
  exit 1
fi

run_one() {
  local m="$1"
  echo "Validating $m"
  "$JAVA" -jar "$JAR" -context reads -manifest "$m" \
    -userName "$ENA_WEBIN_USER" -password "$ENA_WEBIN_PASSWORD" \
    -outputDir "$OUTDIR" -validate
  if [[ "$MODE" == "submit" ]]; then
    echo "Submitting $m"
    "$JAVA" -jar "$JAR" -context reads -manifest "$m" \
      -userName "$ENA_WEBIN_USER" -password "$ENA_WEBIN_PASSWORD" \
      -outputDir "$OUTDIR" -submit
  fi
}

shopt -s nullglob
for m in manifests/r2310_*.manifest.txt; do
  run_one "$m"
done
