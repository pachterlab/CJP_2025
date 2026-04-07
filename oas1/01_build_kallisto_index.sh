#!/usr/bin/env bash
set -euo pipefail

GENCODE_RELEASE="${GENCODE_RELEASE:-49}"
OUT_INDEX="${OUT_INDEX:-gencode.v${GENCODE_RELEASE}.idx}"
SCRATCH_DIR="${SCRATCH_DIR:-}"
TRANSCRIPTS_FA="${TRANSCRIPTS_FA:-}"

if [[ -z "${SCRATCH_DIR}" ]]; then
  SCRATCH_DIR="$(mktemp -d -p . scratch.XXXX)"
  CLEAN_SCRATCH=1
else
  CLEAN_SCRATCH=0
fi

cleanup() {
  if [[ "${CLEAN_SCRATCH}" == "1" ]]; then
    rm -rf "${SCRATCH_DIR}"
  fi
}
trap cleanup EXIT

LOCAL_FA="${TRANSCRIPTS_FA}"
TRANSCRIPTS_FA_TMP="${SCRATCH_DIR}/gencode.v${GENCODE_RELEASE}.transcripts.fa"
TRANSCRIPTS_FA_GZ="${TRANSCRIPTS_FA_TMP}.gz"

BASE_URL="https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_${GENCODE_RELEASE}"
FA_URL="${BASE_URL}/gencode.v${GENCODE_RELEASE}.transcripts.fa.gz"

if [[ -f "${OUT_INDEX}" ]]; then
  echo "Index already exists: ${OUT_INDEX}"
  exit 0
fi

if [[ -n "${LOCAL_FA}" ]]; then
  if [[ ! -f "${LOCAL_FA}" ]]; then
    echo "TRANSCRIPTS_FA not found: ${LOCAL_FA}" >&2
    exit 1
  fi
  kallisto index -i "${OUT_INDEX}" "${LOCAL_FA}"
else
  echo "Downloading ${FA_URL}"
  curl -L --fail "${FA_URL}" -o "${TRANSCRIPTS_FA_GZ}"

  if command -v pigz >/dev/null 2>&1; then
    pigz -dc "${TRANSCRIPTS_FA_GZ}" > "${TRANSCRIPTS_FA_TMP}"
  else
    gzip -dc "${TRANSCRIPTS_FA_GZ}" > "${TRANSCRIPTS_FA_TMP}"
  fi

  kallisto index -i "${OUT_INDEX}" "${TRANSCRIPTS_FA_TMP}"
fi

echo "Built index: ${OUT_INDEX}"
