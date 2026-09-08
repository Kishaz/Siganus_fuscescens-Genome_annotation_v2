#!/usr/bin/env bash
# Submit the pipeline as a Slurm DAG.
#
#   ./scripts/run_all.sh              submit everything
#   DRYRUN=1 ./scripts/run_all.sh     print the plan, submit nothing
#   ./scripts/run_all.sh --from 09    resume from a stage
#
# v2 is not a chain. Helixer needs only the masked genome, so it runs on a GPU
# node while BRAKER3 is still working through RNA-seq and proteins on CPU;
# Liftoff and ncRNA are likewise independent. Submitting a linear chain would
# serialise roughly two days of work for no reason.
#
#   01 ──> 02 repeats ─┬──> 04 rnaseq ──┐
#                      ├──> 05 proteins ┴──> 06 braker3 ──┐
#                      ├──> 07 helixer (GPU) ─────────────┴──> 09 combine
#                      ├──> 08 liftoff ───────────────────────────┐
#                      └──> 11 ncrna ─────────────────┐           │
#   03 fetch (login) ──┘                              │           │
#                                        09 ──> 10 filter <────────┘
#                                             10 ──> 12 functional ──> 13 qc
set -euo pipefail
cd "$(dirname "$0")/.."
source ./config.sh

FROM="00"
while [[ "${1:-}" == --* ]]; do
  case "$1" in
    --from) FROM="$2"; shift 2 ;;
    *) echo "unknown flag: $1" >&2; exit 2 ;;
  esac
done

skip() { [[ "$1" < "$FROM" ]]; }

run_local() {  # run_local <stage> <script>
  if skip "$1"; then echo "[skip] $2"; return; fi
  echo "=== $2 (login node) ==="
  if [ -n "${DRYRUN:-}" ]; then echo "  would run scripts/$2"; return; fi
  bash "scripts/$2"
}

# submit <stage> <script> [dep-jobid ...] -> echoes the new job id
submit() {
  local stage="$1" script="$2"; shift 2
  if skip "$stage"; then echo ""; return; fi
  local deps="" d
  for d in "$@"; do [ -n "$d" ] && deps="${deps:+$deps:}$d"; done
  local args=()
  [ -n "$deps" ] && args+=(--dep "afterok:$deps")
  if [ -n "${DRYRUN:-}" ]; then
    echo "  ${script}${deps:+  after $deps}" >&2
    echo "JOB_$stage"
    return
  fi
  local out id
  out=$(./scripts/submit.sh "${args[@]}" "$script")
  id=$(echo "$out" | grep -oE '[0-9]+$')
  [ -n "$id" ] || { echo "FATAL: no job id from: $out" >&2; exit 1; }
  echo "  $script -> $id${deps:+  (after $deps)}" >&2
  echo "$id"
}

run_local 00 00_setup.sh
run_local 01 01_prepare_genome.sh

echo "=== submitting compute stages ==="
J02=$(submit 02 02_repeats.sbatch)

# Evidence fetching needs the network and so cannot be a job; it runs on the
# login node while RepeatModeler/EarlGrey is queued or grinding away.
run_local 03 03_fetch_evidence.sh

J04=$(submit 04 04_rnaseq.sbatch            "$J02")
J05=$(submit 05 05_protein_evidence.sbatch  "$J02")
J06=$(submit 06 06_braker3.sbatch           "$J04" "$J05")
J07=$(submit 07 07_helixer.sbatch           "$J02")
J08=$(submit 08 08_liftoff.sbatch           "$J02")
J09=$(submit 09 09_combine.sbatch           "$J06" "$J07")
J11=$(submit 11 11_ncrna.sbatch             "$J02")
J10=$(submit 10 10_filter_models.sbatch     "$J09" "$J08")
J12=$(submit 12 12_functional.sbatch        "$J10")
J13=$(submit 13 13_qc.sbatch                "$J12")

cat <<TXT

=== submitted ===
Watch:   squeue -u \$USER
Logs:    $LOGS
Final:   ./scripts/14_package_release.sh    (after job ${J13:-<qc>} succeeds)

Critical path, cpu partition, 48 cores:
  02 repeats (EarlGrey)   3-7 days   <-- dominates everything
  04 rnaseq               4-8 h      }
  05 proteins             2-4 h      } parallel after 02
  07 helixer (GPU)        4-8 h      }
  08 liftoff              2-4 h      }
  11 ncrna                4-8 h      }
  06 braker3              1-2 days
  09 combine              1-2 h
  10 filter               1-2 h
  12 functional           1-2 days   (InterProScan dominates)
  13 qc                   2-4 h

Set REPEAT_METHOD=manual in config.sh to trade EarlGrey's curation for a
faster RepeatModeler2 + RepeatMasker run.
TXT
