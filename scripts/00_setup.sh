#!/usr/bin/env bash
# Stage 00 - provision everything. Login node, needs outbound HTTPS.
# ~1-2 h, ~35 GB (containers dominate).
set -euo pipefail
cd "$(dirname "$0")/.."
source ./config.sh

mkdir -p "$PROJ/bin" "$MAMBA_ROOT_PREFIX" "$SIF_DIR" "$REPO/docs/provenance"

# ---- pick a package manager -------------------------------------------------
# Preference order:
#   1. micromamba already installed
#   2. micromamba staged by mirror_fetch.sh (18 MB) - micro.mamba.pm is blocked
#      by some site proxies, DISCOVERY among them (403 after CONNECT)
#   3. download micromamba, if that host happens to be reachable
#   4. a site conda/mamba
BACKEND=""
if [ -x "$MICROMAMBA_BIN" ]; then
  BACKEND=micromamba
elif [ -x "$STAGED/micromamba" ]; then
  echo "[00] using staged micromamba"
  mkdir -p "$(dirname "$MICROMAMBA_BIN")"
  cp "$STAGED/micromamba" "$MICROMAMBA_BIN" && chmod +x "$MICROMAMBA_BIN"
  BACKEND=micromamba
elif curl -sSfL --max-time 60 https://micro.mamba.pm/api/micromamba/linux-64/latest \
       -o /tmp/mm.$$.tar.bz2 2>/dev/null; then
  echo "[00] downloading micromamba"
  mkdir -p "$(dirname "$MICROMAMBA_BIN")"
  tar -xjf /tmp/mm.$$.tar.bz2 -O bin/micromamba > "$MICROMAMBA_BIN"
  chmod +x "$MICROMAMBA_BIN"; rm -f /tmp/mm.$$.tar.bz2
  BACKEND=micromamba
else
  # Try to surface a site conda via Lmod. On DISCOVERY the anaconda3 module is
  # hidden behind `shared`, so a bare `module load anaconda3` fails.
  if ! command -v conda >/dev/null 2>&1 && command -v module >/dev/null 2>&1; then
    for spec in "shared anaconda3/2023.09" "shared anaconda3" "anaconda3" "anaconda" "miniconda3"; do
      # shellcheck disable=SC1090
      module load $spec >/dev/null 2>&1 && command -v conda >/dev/null 2>&1 && {
        echo "[00] module load $spec"; break; }
    done
  fi
fi

if [ -z "$BACKEND" ]; then
  if command -v mamba >/dev/null 2>&1 || command -v conda >/dev/null 2>&1; then
    echo "[00] micro.mamba.pm unreachable - using the site conda"
    BACKEND=conda
  fi
fi

if [ -z "$BACKEND" ]; then
  cat >&2 <<'TXT'
[00] FATAL: no package manager available.
     micro.mamba.pm is unreachable and no conda is on PATH. Either:
       - module load anaconda   (or whatever your site calls it), or
       - stage the binary:  scripts/mirror_fetch.sh <dir>  then rsync
         <dir>/micromamba into $STAGED/
TXT
  exit 1
fi
echo "[00] backend: $BACKEND"
[ "$BACKEND" = micromamba ] && "$MICROMAMBA_BIN" --version

# ---- environments -----------------------------------------------------------
# Channels come from the yaml files only. --override-channels keeps the
# Anaconda `defaults` channel out of the solve: its terms are a poor fit for a
# public pipeline, and mixing it with conda-forge/bioconda causes solver churn.
SOLVER="$(command -v mamba || command -v conda || true)"
for y in env/*.yaml; do
  n="$(basename "$y" .yaml)"
  envdir="$MAMBA_ROOT_PREFIX/envs/$n"
  if [ -d "$envdir" ]; then
    echo "[00] env '$n' exists - skipping"
    continue
  fi
  echo "[00] creating env '$n'"
  if [ "$BACKEND" = micromamba ]; then
    "$MICROMAMBA_BIN" create -y -r "$MAMBA_ROOT_PREFIX" -f "$y"
  else
    # conda refuses a yaml carrying `name:` when -p is given, so strip it.
    tmpy="$(mktemp)"; grep -v '^name:' "$y" > "$tmpy"
    # `conda env create` is the documented route but 23.x is prone to plugin
    # errors on some site installs; `conda create` with an explicit package
    # list is the blunter, more reliable fallback. Channels come from the
    # yaml/CLI only - this site already has conda-forge + bioconda configured
    # and no `defaults`, so no --override-channels is needed.
    "$SOLVER" env create -q -p "$envdir" -f "$tmpy" \
      || CONDA_NO_PLUGINS=true "$SOLVER" env create -q -p "$envdir" -f "$tmpy" \
      || CONDA_NO_PLUGINS=true "$SOLVER" create -y -q -p "$envdir" \
           -c conda-forge -c bioconda \
           $(awk '/^dependencies:/{f=1;next} f&&/^ *- /{sub(/^ *- /,"");printf "%s ",$0}' "$y") \
      || { echo "[00] FAILED to create env '$n' - see the error above" >&2; rm -f "$tmpy"; continue; }
    rm -f "$tmpy"
  fi
  [ -d "$envdir" ] && echo "[00]   -> $envdir"
done

# ---- freeze exact solved versions ------------------------------------------
# `>=` constraints in the yaml files are what we asked for; this records what we
# actually got, so the run is reproducible even after upstream moves on.
echo "[00] freezing resolved environments -> docs/provenance/"
for y in env/*.yaml; do
  n="$(basename "$y" .yaml)"
  "$MICROMAMBA_BIN" env export -r "$MAMBA_ROOT_PREFIX" -n "$n" --explicit \
    > "docs/provenance/${n}.lock.txt" 2>/dev/null \
  || "$MICROMAMBA_BIN" list -r "$MAMBA_ROOT_PREFIX" -n "$n" \
    > "docs/provenance/${n}.versions.txt" 2>/dev/null || true
done

# ---- containers -------------------------------------------------------------
if ! command -v "$CONTAINER_CMD" >/dev/null; then
  echo "[00] FATAL: $CONTAINER_CMD not found; BRAKER3 and Helixer both need it" >&2
  exit 1
fi

pull() {  # pull <image-uri> <sif-name>
  local uri="$1" sif="$SIF_DIR/$2"
  if [ -s "$sif" ]; then echo "[00] have $2"; return; fi
  echo "[00] building $2 from $uri"
  "$CONTAINER_CMD" build "$sif" "$uri"
}

pull "$BRAKER_IMG" braker3.sif
if [ "${USE_HELIXER:-1}" = "1" ]; then
  pull "$HELIXER_IMG" helixer.sif
fi

# Record what the images actually are - the tag alone is not a fingerprint.
{
  echo -e "image\turi\tsif_md5\tbuilt"
  for pair in "braker3.sif|$BRAKER_IMG" "helixer.sif|$HELIXER_IMG"; do
    f="${pair%%|*}"; u="${pair##*|}"
    [ -s "$SIF_DIR/$f" ] || continue
    printf '%s\t%s\t%s\t%s\n' "$f" "$u" \
      "$(md5sum "$SIF_DIR/$f" | awk '{print $1}')" "$(date -Is)"
  done
} > docs/provenance/containers.tsv
cat docs/provenance/containers.tsv

echo
echo "[00] done. Next: ./scripts/01_prepare_genome.sh"
