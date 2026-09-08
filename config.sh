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
export PROJ="${PROJ:-/projects/yates_lab_hpc/sam/smwambu1/aigo_genome}"
export WORK="$PROJ/work"          # intermediates (large, purgeable)
export DB="$PROJ/db"              # reference databases (large, reusable)
export REL="$PROJ/release"        # final deliverables
export LOGS="$PROJ/logs"
export REPO="${REPO:-$PROJ/Siganus_fuscescens-Genome_annotation_v2}"
export LIB="$REPO/lib"
export STAGED="$DB/staged"        # firewall-blocked downloads, rsynced in

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
mm() {  # run a command inside a micromamba env:  mm <env> <cmd...>
  "$MICROMAMBA_BIN" run -r "$MAMBA_ROOT_PREFIX" -n "$1" "${@:2}"
}
need() { for f in "$@"; do [ -s "$f" ] || { echo "MISSING INPUT: $f" >&2; exit 1; }; done; }
done_stamp() { mkdir -p "$WORK/.stamps"; touch "$WORK/.stamps/$1"; }
have_stamp() { [ -f "$WORK/.stamps/$1" ]; }
export -f mm need done_stamp have_stamp

mkdir -p "$WORK" "$DB" "$REL" "$LOGS" "$SIF_DIR" "$WORK/tmp" "$STAGED" 2>/dev/null || true
