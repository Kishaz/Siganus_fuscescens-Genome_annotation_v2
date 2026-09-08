#!/usr/bin/env bash
# ============================================================================
# Siganus fuscescens (Sfus_1.0) annotation v2 - site configuration
#
# Edit ONLY this file to move the pipeline to a different cluster. No
# #SBATCH --account or --partition is hard-coded anywhere; scripts/submit.sh
# supplies them from here.
# ============================================================================

# --- Slurm ------------------------------------------------------------------
export SLURM_ACCOUNT="${SLURM_ACCOUNT:-yates_lab_hpc}"    # or yateslab_genomics
export PART_CPU="${PART_CPU:-cpu}"          # 143 nodes, 48+ cores, 488G+, 10 d
export PART_BIGMEM="${PART_BIGMEM:-cpu}"    # cpu nodes already carry 488G+
export PART_GPU="${PART_GPU:-gpuh100}"      # 4x H100, 128 cores - Helixer
export MAIL_USER="${MAIL_USER:-}"           # blank = no mail

# --- Paths ------------------------------------------------------------------
# Two tiers, deliberately:
#   SCRATCH  big, fast, PURGEABLE - intermediates and reference databases
#   PROJ     durable, quota'd     - the repo, logs, checkpoints, the release
# Nothing irreplaceable is allowed to live only on scratch. Stage 02 alone is
# 3-7 days of compute, so its outputs are checkpointed to PROJ as they appear.
export PROJ="${PROJ:-/projects/yates_lab_hpc/sam/smwambu1/aigo_genome}"

# SCRATCH is auto-detected rather than assumed. On DISCOVERY as of 2026-09,
# /scratch is 100% full cluster-wide (mkdir returns ENOSPC), so the pipeline
# must not simply trust that it exists. Set SCRATCH explicitly to override.
#
# Needs ~300 GB: ~150 GB databases, ~150 GB intermediates.
: "${SCRATCH_MIN_GB:=350}"
if [ -z "${SCRATCH:-}" ]; then
  for cand in "/scratch/$USER/aigo_genome" "/lscratch/$USER/aigo_genome"; do
    mkdir -p "$cand" 2>/dev/null || continue
    # Never site 300 GB of intermediates on a RAM-backed filesystem: on this
    # cluster the compute-node root is tmpfs, and filling it exhausts node
    # memory rather than disk.
    fstype=$(df -PT "$cand" 2>/dev/null | tail -1 | awk '{print $2}')
    case "$fstype" in
      tmpfs|ramfs|devtmpfs)
        echo "NOTE: skipping $cand - it is $fstype (RAM-backed)." >&2
        rmdir "$cand" 2>/dev/null || true; continue ;;
    esac
    avail_gb=$(df -BG --output=avail "$cand" 2>/dev/null | tail -1 | tr -dc '0-9')
    if [ -n "$avail_gb" ] && [ "$avail_gb" -ge "$SCRATCH_MIN_GB" ]; then
      SCRATCH="$cand"; break
    fi
    rmdir "$cand" 2>/dev/null || true
  done
fi
if [ -z "${SCRATCH:-}" ]; then
  SCRATCH="$PROJ/scratch"
  mkdir -p "$SCRATCH" 2>/dev/null || true
  if [ "${SCRATCH_WARNED:-0}" != "1" ]; then
    echo "NOTE: no scratch filesystem with >=${SCRATCH_MIN_GB} GB free." >&2
    echo "      Using $SCRATCH on /projects, which counts against the group" >&2
    echo "      quota. Set SCRATCH=... in config.sh once /scratch has room." >&2
    export SCRATCH_WARNED=1
  fi
fi
export SCRATCH

export WORK="$SCRATCH/work"       # intermediates (large, regenerable)
export DB="$SCRATCH/db"           # reference databases (large, re-downloadable)
export STAGED="$DB/staged"        # firewall-blocked downloads, rsynced in

export REL="$PROJ/release"        # final deliverables - durable
export LOGS="$PROJ/logs"          # durable, so a failure is diagnosable later
export CKPT="$PROJ/checkpoints"   # expensive intermediates, survive a purge
export REPO="${REPO:-$PROJ/Siganus_fuscescens-Genome_annotation_v2}"
export LIB="$REPO/lib"

# Node-local fast scratch; falls back to $WORK/tmp.
export FASTTMP="${SLURM_TMPDIR:-${TMPDIR:-$WORK/tmp}}"

# --- Software provisioning --------------------------------------------------
export MAMBA_ROOT_PREFIX="${MAMBA_ROOT_PREFIX:-$PROJ/micromamba}"
export MICROMAMBA_BIN="${MICROMAMBA_BIN:-$PROJ/bin/micromamba}"
export CONTAINER_CMD="${CONTAINER_CMD:-apptainer}"
export SIF_DIR="$PROJ/containers"

# Containers are pinned to explicit versions, never :latest - a re-run six
# months from now must reproduce this run, and BRAKER/Helixer both move fast.
export BRAKER_IMG="docker://teambraker/braker3:v3.1.1"
export HELIXER_IMG="docker://gglyptodon/helixer-docker:helixer_v0.3.7_cuda_12.2.2-cudnn8"

# --- Assembly ---------------------------------------------------------------
export ASM_ACC="GCA_042847805.1"
export ASM_NAME="Sfus_1.0"
export SPECIES="Siganus_fuscescens"
export TAXID="225757"
export ANNOT_VERSION="v2.0"
export GENE_PREFIX="SFUS"          # stable IDs: SFUS_G000001 / SFUS_G000001.1

# --- QC lineages ------------------------------------------------------------
# odb10 is the GATED lineage: the 2024 paper reported odb10 (98.5%, genome
# mode), so only an odb10 number is comparable. odb12 is reported alongside.
export BUSCO_LINEAGE="actinopterygii_odb10"
export BUSCO_LINEAGES="actinopterygii_odb10 actinopterygii_odb12"

# --- Evidence ---------------------------------------------------------------
# Conspecific RNA-seq only, by project decision. The 80 public S. canaliculatus
# RNA-seq runs are deliberately excluded to keep provenance unambiguous; the
# congener contributes through Liftoff and the protein panel instead.
export RNASEQ_RUNS="SRR36642802"

export ORTHODB_URL="https://bioinf.uni-greifswald.de/bioinf/partitioned_odb12/Vertebrata.fa.gz"
# BRAKER wants the smallest OrthoDB clade containing the target. v12 ships no
# Actinopterygii partition and Vertebrata is 19.4 M proteins (68% tetrapod).
# Headers carry the NCBI taxid, so we cut the clade ourselves: 6.17 M proteins
# from 237 ray-finned species, a 3.1x reduction. Set 0 to use all of Vertebrata.
export ORTHODB_CLADE_TAXID="${ORTHODB_CLADE_TAXID:-7898}"    # Actinopterygii
export TAXDUMP_URL="https://ftp.ncbi.nlm.nih.gov/pub/taxonomy/taxdump.tar.gz"

# Congener with a chromosome-level annotated genome - used by Liftoff as an
# independent evidence track and as a synteny sanity check.
export LIFTOFF_REF_ACC="GCA_053572195.1"     # Siganus canaliculatus, QU_Sigcan_1.0

export DFAM_BASE="https://www.dfam.org/releases/current/families/FamDB"
export DFAM_FILES="dfam40.0.h5.gz dfam40.curated.consensus.0.h5.gz"
export DFAM_RELEASE="4.0"

# --- Repeat strategy --------------------------------------------------------
# "earlgrey"  EarlGrey 7.3.1 - RepeatModeler2 plus BLAST-extend-extend curation,
#             TE classification and a landscape plot. Current best practice and
#             what v2 uses by default.
# "manual"    the v1 route (RepeatModeler2 -> TEsorter -> 2-pass RepeatMasker),
#             kept because it is faster and needs no extra dependencies.
export REPEAT_METHOD="${REPEAT_METHOD:-earlgrey}"

# --- Gene-model combination -------------------------------------------------
# BRAKER3 is the evidence-based backbone. Helixer is an independent
# deep-learning prediction that uses no evidence at all, so where the two agree
# the model is supported by genuinely orthogonal reasoning. TSEBRA merges them.
export USE_HELIXER="${USE_HELIXER:-1}"
export HELIXER_LINEAGE="${HELIXER_LINEAGE:-vertebrate}"

# --- Resources --------------------------------------------------------------
export THREADS_BIG="${THREADS_BIG:-48}"
export MEM_BIG="${MEM_BIG:-400G}"
export TIME_LONG="${TIME_LONG:-7-00:00:00}"

# --- Helpers ----------------------------------------------------------------
# Keep conda's package cache off /home. The site default is ~/.conda/pkgs, and
# on this cluster /home is a 98 GB volume with ~18 GB free - eight environments
# of bioinformatics packages will exhaust it and the failure mode is confusing.
export CONDA_PKGS_DIRS="${CONDA_PKGS_DIRS:-$SCRATCH/conda_pkgs}"
export MAMBA_PKGS_DIRS="$CONDA_PKGS_DIRS"
mkdir -p "$CONDA_PKGS_DIRS" 2>/dev/null || true

# conda_sh - locate the profile script that makes `conda activate` work in a
# non-interactive shell. `conda activate` fails without it, which is why
# `source activate` is often the only thing that works on a site install.
conda_sh() {
  local base
  base="$(conda info --base 2>/dev/null)" || return 1
  [ -s "$base/etc/profile.d/conda.sh" ] && { echo "$base/etc/profile.d/conda.sh"; return 0; }
  return 1
}

# mm <env> <cmd...> - run a command inside a pipeline environment.
#
# Backend-agnostic on purpose. micro.mamba.pm is blocked by some site proxies
# (DISCOVERY returns 403 after CONNECT), so micromamba cannot be assumed
# installable. Both back ends place envs at $MAMBA_ROOT_PREFIX/envs/<name>, so
# the layout is identical either way.
mm() {
  local env="$1"; shift
  local envdir="$MAMBA_ROOT_PREFIX/envs/$env"

  if [ -x "$MICROMAMBA_BIN" ]; then
    "$MICROMAMBA_BIN" run -r "$MAMBA_ROOT_PREFIX" -n "$env" "$@"
    return $?
  fi

  if ! command -v conda >/dev/null 2>&1; then
    echo "FATAL: no backend for env '$env'. Stage micromamba into \$STAGED," >&2
    echo "       or: module load shared anaconda3/2023.09" >&2
    return 1
  fi

  # `conda run` is cleanest but is fragile on older conda (23.x hits plugin and
  # locking bugs under load). Fall back to activating in a subshell, which is
  # what `source activate` does and is the reliable path on this site install.
  if conda run -p "$envdir" --no-capture-output true >/dev/null 2>&1; then
    conda run -p "$envdir" --no-capture-output "$@"
  else
    local csh; csh="$(conda_sh || true)"
    if [ -n "$csh" ]; then
      ( set +u; . "$csh"; conda activate "$envdir" || return 1; exec "$@" )
    else
      ( set +u; . activate "$envdir" 2>/dev/null || return 1; exec "$@" )
    fi
  fi
}
need() { for f in "$@"; do [ -s "$f" ] || { echo "MISSING INPUT: $f" >&2; exit 1; }; done; }

# Stamps live beside the data they describe. If scratch is purged the stamps go
# with it, so a rerun correctly redoes the work rather than trusting a stale
# marker for files that no longer exist.
done_stamp() { mkdir -p "$WORK/.stamps"; touch "$WORK/.stamps/$1"; }
have_stamp() { [ -f "$WORK/.stamps/$1" ]; }

# checkpoint <file> [<file>...] - copy an expensive artefact to durable storage.
# Used for anything that costs more to regenerate than to store: the curated
# repeat library, the softmasked genome, raw predictions.
checkpoint() {
  mkdir -p "$CKPT"
  local f
  for f in "$@"; do
    [ -s "$f" ] || continue
    if [ ! -e "$CKPT/$(basename "$f")" ] || [ "$f" -nt "$CKPT/$(basename "$f")" ]; then
      cp -a "$f" "$CKPT/" && echo "  checkpointed $(basename "$f") -> $CKPT"
    fi
  done
}

# restore_checkpoint <dest-dir> <file> - bring one back after a scratch purge.
restore_checkpoint() {
  local dest="$1" f="$2"
  if [ ! -s "$dest/$f" ] && [ -s "$CKPT/$f" ]; then
    mkdir -p "$dest" && cp -a "$CKPT/$f" "$dest/" && echo "  restored $f from checkpoint"
  fi
}
export -f mm need done_stamp have_stamp checkpoint restore_checkpoint

mkdir -p "$WORK" "$DB" "$REL" "$LOGS" "$CKPT" "$SIF_DIR" "$WORK/tmp" "$STAGED" 2>/dev/null || true
