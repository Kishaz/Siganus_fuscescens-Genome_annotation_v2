#!/usr/bin/env bash
# Portable submit wrapper: keeps account/partition out of the .sbatch files so
# the pipeline moves between clusters by editing config.sh only.
#
#   ./scripts/submit.sh 02_repeatmodeler.sbatch
#   ./scripts/submit.sh --dep afterok:12345 03_repeatmasker.sbatch
#   DRYRUN=1 ./scripts/submit.sh 07_braker3.sbatch     # print, do not submit
set -euo pipefail
cd "$(dirname "$0")/.."
source ./config.sh

DEP=""
while [[ "${1:-}" == --* ]]; do
  case "$1" in
    --dep) DEP="--dependency=$2"; shift 2 ;;
    --part) PART_CPU="$2"; shift 2 ;;
    *) echo "unknown flag $1" >&2; exit 2 ;;
  esac
done

JOB="${1:?usage: submit.sh [--dep afterok:JOBID] <script.sbatch>}"
[ -f "scripts/$JOB" ] || JOB="$(basename "$JOB")"
[ -f "scripts/$JOB" ] || { echo "no such script: scripts/$JOB" >&2; exit 1; }

# Per-job partition override: a script may declare  #PIPE-PARTITION: gpu
P="$(grep -m1 '^#PIPE-PARTITION:' "scripts/$JOB" | awk '{print $2}' || true)"
case "$P" in
  gpu)    PART="$PART_GPU" ;;
  bigmem) PART="$PART_BIGMEM" ;;
  *)      PART="$PART_CPU" ;;
esac

MAILOPT=()
[ -n "${MAIL_USER:-}" ] && MAILOPT=(--mail-user="$MAIL_USER" --mail-type=END,FAIL)

mkdir -p "$LOGS"
CMD=(sbatch -A "$SLURM_ACCOUNT" -p "$PART"
     -o "$LOGS/%x.%j.out" -e "$LOGS/%x.%j.err"
     --export=ALL,REPO="$REPO",PROJ="$PROJ"
     "${MAILOPT[@]}" ${DEP:+$DEP} "scripts/$JOB")

if [ -n "${DRYRUN:-}" ]; then printf '%q ' "${CMD[@]}"; echo; exit 0; fi
"${CMD[@]}"
