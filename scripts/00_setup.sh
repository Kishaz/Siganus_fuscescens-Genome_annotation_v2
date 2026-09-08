#!/usr/bin/env bash
# Stage 00 - provision everything. Login node, needs outbound HTTPS.
# ~1-2 h, ~35 GB (containers dominate).
set -euo pipefail
cd "$(dirname "$0")/.."
source ./config.sh

mkdir -p "$PROJ/bin" "$MAMBA_ROOT_PREFIX" "$SIF_DIR" "$REPO/docs/provenance"

# ---- micromamba -------------------------------------------------------------
# A static binary in project space: no root, and independent of the site
# anaconda install (whose `defaults` channel carries licensing conditions we
# have no reason to accept for a public pipeline).
if [ ! -x "$MICROMAMBA_BIN" ]; then
  echo "[00] installing micromamba -> $MICROMAMBA_BIN"
  curl -sSL --retry 5 https://micro.mamba.pm/api/micromamba/linux-64/latest \
    | tar -xj -C "$PROJ" bin/micromamba
  chmod +x "$MICROMAMBA_BIN"
fi
"$MICROMAMBA_BIN" --version

# ---- environments -----------------------------------------------------------
for y in env/*.yaml; do
  n="$(basename "$y" .yaml)"
  if "$MICROMAMBA_BIN" env list -r "$MAMBA_ROOT_PREFIX" 2>/dev/null \
       | awk '{print $1}' | grep -qx "$n"; then
    echo "[00] env '$n' exists - skipping"
  else
    echo "[00] creating env '$n'"
    "$MICROMAMBA_BIN" create -y -r "$MAMBA_ROOT_PREFIX" -f "$y"
  fi
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
